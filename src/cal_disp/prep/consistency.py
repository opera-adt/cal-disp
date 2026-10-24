"""Geometric consistency of the calibration inputs.

The workflow combines the DISP displacement, the static LOS vectors and the
DEM by array index, so they must be on the same grid: equal shape, transform
and CRS, and of the same OPERA frame. Static layers of another frame with the
same shape would otherwise silently project the GNSS field onto the wrong
line of sight. :func:`check_input_geometry` verifies this once, before any
processing, and raises a ``ValueError`` naming the offending file and values.

LOS pixels without data (all-zero vectors, or the raster's nodata value) are
set to NaN by :func:`mask_los_nodata`, so that they cannot be used as a valid
line of sight.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import rasterio
import xarray as xr
from rasterio.crs import CRS
from rasterio.transform import Affine

from cal_disp.product import DispProduct, StaticLayer

logger = logging.getLogger(__name__)

#: Largest allowed difference of each transform term, in pixels.
PIXEL_TOLERANCE: float = 1e-6


@dataclass(frozen=True)
class GridGeometry:
    """Grid definition of one input raster."""

    label: str
    path: Path
    shape: tuple[int, int]
    transform: Affine
    crs: CRS | None
    frame_id: int | None

    def describe(self) -> str:
        """One-line summary of the compared values."""
        t = self.transform
        epsg = self.crs.to_epsg() if self.crs else None
        crs = f"EPSG:{epsg}" if epsg else (self.crs.to_string() if self.crs else "none")
        frame = f"F{self.frame_id:05d}" if self.frame_id is not None else "unknown"
        return (
            f"{self.label}: frame {frame}, shape {self.shape}, {crs},"
            f" origin ({t.c:.3f}, {t.f:.3f}), pixel ({t.a:g}, {t.e:g})"
        )


def disp_geometry(disp: DispProduct | Path | str) -> GridGeometry:
    """Grid of a DISP product.

    The transform is the ``GeoTransform`` of ``spatial_ref`` when present,
    otherwise it is derived from the (pixel-centre) ``x``/``y`` coordinates.
    """
    if not isinstance(disp, DispProduct):
        disp = DispProduct.from_path(disp)
    with disp.open_dataset() as ds:
        shape = (ds.sizes["y"], ds.sizes["x"])
        attrs = ds["spatial_ref"].attrs if "spatial_ref" in ds else {}
        crs_wkt = attrs.get("crs_wkt")
        geotransform = attrs.get("GeoTransform")
        if geotransform is not None:
            c, a, b, f, d, e = (float(v) for v in str(geotransform).split())
            transform = Affine(a, b, c, d, e, f)
        else:
            x, y = ds["x"].values, ds["y"].values
            if len(x) < 2 or len(y) < 2:
                raise ValueError(
                    f"DISP file {disp.path} has no GeoTransform and fewer than 2"
                    " pixels per axis: cannot determine its grid"
                )
            dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
            transform = Affine(
                dx, 0.0, float(x[0]) - dx / 2, 0.0, dy, float(y[0]) - dy / 2
            )
    return GridGeometry(
        label="DISP",
        path=disp.path,
        shape=shape,
        transform=transform,
        crs=CRS.from_wkt(crs_wkt) if crs_wkt else None,
        frame_id=disp.frame_id,
    )


def raster_geometry(path: Path | str, label: str) -> GridGeometry:
    """Grid of a GeoTIFF; the frame id is parsed from an OPERA static name."""
    path = Path(path)
    with rasterio.open(path) as src:
        shape, transform, crs = (src.height, src.width), src.transform, src.crs
    try:
        frame_id: int | None = StaticLayer.from_path(path).frame_id
    except ValueError:
        frame_id = None
    return GridGeometry(label, path, shape, transform, crs, frame_id)


def _grid_mismatches(
    ref: GridGeometry, other: GridGeometry, pixel_tolerance: float
) -> list[str]:
    """Differences between two grids, as readable statements."""
    problems: list[str] = []
    if other.frame_id is not None and other.frame_id != ref.frame_id:
        problems.append(f"frame id {other.frame_id} != {ref.frame_id}")
    if other.shape != ref.shape:
        problems.append(f"shape {other.shape} != {ref.shape}")
    if ref.crs is None or other.crs is None or other.crs != ref.crs:
        problems.append(f"CRS {_crs_name(other.crs)} != {_crs_name(ref.crs)}")

    # Each transform term within `pixel_tolerance` pixels
    tol_x = abs(ref.transform.a) * pixel_tolerance
    tol_y = abs(ref.transform.e) * pixel_tolerance
    terms = ("a", "b", "c", "d", "e", "f")
    tols = (tol_x, tol_x, tol_x, tol_y, tol_y, tol_y)
    bad = [
        name
        for name, tol in zip(terms, tols)
        if abs(getattr(other.transform, name) - getattr(ref.transform, name)) > tol
    ]
    if bad:
        o, r = other.transform, ref.transform
        shift_x = (o.c - r.c) / r.a if r.a else float("nan")
        shift_y = (o.f - r.f) / r.e if r.e else float("nan")
        problems.append(
            f"transform (origin ({o.c}, {o.f}), pixel ({o.a}, {o.e})) != (origin"
            f" ({r.c}, {r.f}), pixel ({r.a}, {r.e})); origin shifted by"
            f" ({shift_x:g}, {shift_y:g}) pixels"
        )
    return problems


def _crs_name(crs: CRS | None) -> str:
    if crs is None:
        return "none"
    epsg = crs.to_epsg()
    return f"EPSG:{epsg}" if epsg else crs.to_string()


def check_tropo_coverage(
    disp: DispProduct | Path | str, tropo_files: Sequence[Path | str]
) -> None:
    """Raise ``ValueError`` unless each TROPO file covers the DISP frame."""
    if not tropo_files:
        return
    if not isinstance(disp, DispProduct):
        disp = DispProduct.from_path(disp)
    b = disp.get_bounds_wgs84()
    for tropo_file in tropo_files:
        with xr.open_dataset(tropo_file, engine="h5netcdf") as ds:
            if "latitude" not in ds.coords or "longitude" not in ds.coords:
                raise ValueError(
                    f"TROPO file {tropo_file} has no latitude/longitude coordinates"
                )
            lat, lon = ds["latitude"].values, ds["longitude"].values
        south, north = float(lat.min()), float(lat.max())
        west, east = float(lon.min()), float(lon.max())
        if (
            south > b["south"]
            or north < b["north"]
            or west > b["west"]
            or east < b["east"]
        ):
            raise ValueError(
                f"TROPO file {tropo_file} does not cover DISP file"
                f" {disp.path.name}: TROPO extent (W {west:.3f}, S {south:.3f},"
                f" E {east:.3f}, N {north:.3f}) vs frame (W {b['west']:.3f},"
                f" S {b['south']:.3f}, E {b['east']:.3f}, N {b['north']:.3f})"
            )


def check_input_geometry(
    disp: DispProduct | Path | str,
    los_file: Path | str,
    dem_file: Path | str | None = None,
    tropo_files: Sequence[Path | str] = (),
    pixel_tolerance: float = PIXEL_TOLERANCE,
) -> dict[str, GridGeometry]:
    """Verify that the static layers are on the grid of the DISP product.

    Parameters
    ----------
    disp : DispProduct, Path or str
        DISP product defining the reference grid.
    los_file : Path or str
        LOS GeoTIFF (``..._line_of_sight_enu.tif``).
    dem_file : Path or str, optional
        DEM GeoTIFF (``..._dem.tif``).
    tropo_files : sequence of Path or str, optional
        TROPO-ZENITH files; each must cover the DISP frame (they are on a
        geographic grid of their own and are interpolated, not indexed).
    pixel_tolerance : float, optional
        Largest difference of each transform term, in pixels. Default 1e-6.

    Returns
    -------
    dict[str, GridGeometry]
        The compared geometries, keyed by ``"DISP"``, ``"LOS"``, ``"DEM"``.

    Raises
    ------
    ValueError
        If a static layer differs from the DISP grid in shape, transform or
        CRS, belongs to another frame (frame id parsed from OPERA file
        names), or a TROPO file does not cover the frame. The message names
        the file and the differing values.

    """
    ref = disp_geometry(disp)
    if ref.crs is None:
        raise ValueError(f"DISP file {ref.path} has no CRS (spatial_ref crs_wkt)")
    geometries = {"DISP": ref}
    logger.info("Input geometry %s", ref.describe())

    for label, path in (("LOS", los_file), ("DEM", dem_file)):
        if path is None:
            continue
        geom = raster_geometry(path, label)
        problems = _grid_mismatches(ref, geom, pixel_tolerance)
        if problems:
            raise ValueError(
                f"{label} file {geom.path} is not on the grid of DISP file"
                f" {ref.path.name}: "
                + "; ".join(problems)
            )
        if geom.frame_id is None:
            logger.warning(
                "%s file %s is not an OPERA DISP-S1-STATIC name: frame id not checked",
                label,
                geom.path.name,
            )
        geometries[label] = geom
        logger.info("Input geometry %s", geom.describe())

    check_tropo_coverage(disp, tropo_files)
    return geometries


def mask_los_nodata(
    los_east: np.ndarray,
    los_north: np.ndarray,
    los_up: np.ndarray,
    nodata: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Set LOS pixels without data to NaN in all three components.

    A pixel has no data when the vector is all zero (a unit vector never is;
    OPERA static layers use 0 as nodata) or when any component equals the
    raster's `nodata` value.
    """
    invalid = (los_east == 0) & (los_north == 0) & (los_up == 0)
    if nodata is not None and not np.isnan(nodata):
        invalid |= (los_east == nodata) | (los_north == nodata) | (los_up == nodata)
    n_invalid = int(invalid.sum())
    if n_invalid:
        logger.info(
            "LOS: %d of %d pixels without data set to NaN", n_invalid, invalid.size
        )
        los_east, los_north, los_up = (
            np.where(invalid, np.nan, band).astype(band.dtype, copy=False)
            for band in (los_east, los_north, los_up)
        )
    return los_east, los_north, los_up
