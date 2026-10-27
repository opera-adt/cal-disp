"""Utility functions for calibration products."""

from datetime import datetime

import numpy as np
import xarray as xr
from rasterio.transform import Affine

from .._disp import DispProduct


def build_filename(
    disp_product: DispProduct,
    sensor: str,
    version: str,
    production_date: datetime,
) -> str:
    """Build OPERA-compliant filename for calibration product."""
    return (
        f"OPERA_L4_DISP-CAL-{sensor}_"
        f"{disp_product.mode}_"
        f"F{disp_product.frame_id:05d}_"
        f"{disp_product.polarization}_"
        f"{disp_product.reference_date:%Y%m%dT%H%M%S}Z_"
        f"{disp_product.secondary_date:%Y%m%dT%H%M%S}Z_"
        f"v{version}_"
        f"{production_date:%Y%m%dT%H%M%S}Z.nc"
    )


def get_transform(ds: xr.Dataset) -> Affine:
    """Extract affine transform from dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with spatial_ref containing GeoTransform.

    Returns
    -------
    Affine
        Affine transformation.

    Raises
    ------
    ValueError
        If GeoTransform not found.

    """
    gt = ds.spatial_ref.attrs.get("GeoTransform")
    if gt is None:
        raise ValueError("No GeoTransform found in spatial_ref")

    vals = [float(x) for x in gt.split()]
    return Affine(vals[1], vals[2], vals[0], vals[4], vals[5], vals[3])


def get_crs(ds: xr.Dataset) -> str:
    """Extract CRS from dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with spatial_ref containing crs_wkt.

    Returns
    -------
    str
        CRS as WKT string.

    Raises
    ------
    ValueError
        If crs_wkt not found.

    """
    crs_wkt = ds.spatial_ref.attrs.get("crs_wkt")
    if crs_wkt is None:
        raise ValueError("No crs_wkt found in spatial_ref")
    return crs_wkt


def compute_transform_from_coords(x: np.ndarray, y: np.ndarray) -> Affine:
    """Affine transform (upper-left pixel *edge*) from pixel-centre coordinates."""
    dx = float(x[1] - x[0])
    dy = float(y[1] - y[0])
    return Affine.translation(float(x[0]) - dx / 2, float(y[0]) - dy / 2) * (
        Affine.scale(dx, dy)
    )


def geotransform_from_coords(x: np.ndarray, y: np.ndarray) -> str:
    """GDAL ``GeoTransform`` string for a grid of pixel-centre coordinates.

    Matches the DISP-S1 convention (upper-left pixel edge as origin), e.g.
    ``"71970.0 30.0 0.0 3385920.0 0.0 -30.0"``.
    """
    t = compute_transform_from_coords(x, y)
    return f"{t.c} {t.a} {t.b} {t.f} {t.d} {t.e}"


def make_spatial_ref(
    spatial_ref: xr.DataArray, x: np.ndarray, y: np.ndarray
) -> xr.DataArray:
    """``spatial_ref`` variable for the grid (`y`, `x`).

    The CRS attributes of the input (``crs_wkt`` and the CF grid-mapping
    attributes) are kept verbatim; ``GeoTransform`` is set for this grid so
    the variable stays georeferenced when written for a coarser grid too.
    If the input carries only ``crs_wkt``, the CF grid-mapping attributes
    are derived from it.
    """
    attrs = {
        k: v.item() if isinstance(v, np.ndarray) and v.size == 1 else v
        for k, v in spatial_ref.attrs.items()
    }
    if "grid_mapping_name" not in attrs and attrs.get("crs_wkt"):
        import pyproj

        attrs = {**pyproj.CRS.from_wkt(attrs["crs_wkt"]).to_cf(), **attrs}
    attrs["GeoTransform"] = geotransform_from_coords(x, y)
    attrs.setdefault("units", "unitless")
    attrs.setdefault(
        "long_name", "Dummy variable with geo-referencing metadata in attributes"
    )
    return xr.DataArray(np.int64(0), attrs=attrs, name="spatial_ref")


def grid_bounds(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    """Outer pixel edges ``(west, south, east, north)`` of a pixel-centre grid."""
    t = compute_transform_from_coords(x, y)
    west, north = t.c, t.f
    east = west + t.a * len(x)
    south = north + t.e * len(y)
    return (min(west, east), min(south, north), max(west, east), max(south, north))


def bounding_polygon_wkt(
    x: np.ndarray, y: np.ndarray, crs_wkt: str, points_per_edge: int = 10
) -> str:
    """WKT ``POLYGON`` (longitude latitude, WGS 84) of the grid's outer edges.

    Each edge is densified with `points_per_edge` vertices before
    reprojection so the polygon follows the curved UTM boundary.
    """
    from rasterio.warp import transform

    west, south, east, north = grid_bounds(x, y)
    n = max(points_per_edge, 2)
    xs = np.concatenate(
        [
            np.linspace(west, east, n),  # north edge, W -> E
            np.full(n, east),  # east edge, N -> S
            np.linspace(east, west, n),  # south edge, E -> W
            np.full(n, west),  # west edge, S -> N
        ]
    )
    ys = np.concatenate(
        [
            np.full(n, north),
            np.linspace(north, south, n),
            np.full(n, south),
            np.linspace(south, north, n),
        ]
    )
    lons, lats = transform(crs_wkt, "EPSG:4326", xs.tolist(), ys.tolist())
    ring = [f"{lon:.6f} {lat:.6f}" for lon, lat in zip(lons, lats)]
    ring.append(ring[0])
    return "POLYGON((" + ", ".join(ring) + "))"


RASTER_CHUNK = 256


def build_encoding(
    ds: xr.Dataset, compression: bool = True, chunk: int = RASTER_CHUNK
) -> dict[str, dict]:
    """NetCDF encoding for every variable of `ds`.

    Raster variables (2-D ``(y, x)`` or 3-D ``(time, y, x)``) are chunked
    ``(chunk, chunk)`` / ``(1, chunk, chunk)`` and, with `compression`,
    written with gzip level 4 plus the shuffle filter, as the DISP-S1 input.
    Floating-point rasters keep ``float32`` and ``_FillValue = NaN``.
    Coordinate variables get no ``_FillValue`` (CF forbids it on them) and
    no compression; scalar and string variables are left as they are.
    """
    encoding: dict[str, dict] = {}
    for name_, var in ds.data_vars.items():
        name = str(name_)
        if var.ndim < 2:
            continue
        enc: dict = {
            "chunksizes": tuple(
                min(chunk, size) if dim in ("y", "x") else 1
                for dim, size in zip(var.dims, var.shape)
            )
        }
        if compression:
            enc.update(zlib=True, complevel=4, shuffle=True)
        else:
            enc.update(zlib=False, shuffle=False)
        if np.issubdtype(var.dtype, np.floating):
            enc.update(dtype="float32", _FillValue=np.float32(np.nan))
        encoding[name] = enc
    for coord in ds.coords:
        if ds[coord].ndim == 1:
            encoding[str(coord)] = {"_FillValue": None}
    return encoding


def compute_stats(data: np.ndarray) -> dict[str, float] | None:
    """Compute summary statistics for array, ignoring NaN values."""
    valid_data = data[~np.isnan(data)]

    if len(valid_data) == 0:
        return None

    return {
        "mean": float(np.mean(valid_data)),
        "std": float(np.std(valid_data)),
        "min": float(np.min(valid_data)),
        "max": float(np.max(valid_data)),
    }
