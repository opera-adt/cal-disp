"""Calibrate one OPERA DISP-S1 product against GNSS, using Venti.

Load the product, compute the GNSS LOS field from pre-staged UNR files,
estimate the calibration surface (``venti.estimate_calibration_surface``)
and write a ``CalProduct``.
"""

from __future__ import annotations

import functools
import importlib.metadata
import inspect
import json
import logging
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np
import rasterio
import xarray as xr
from rasterio.windows import Window

from cal_disp._version import __version__
from cal_disp.config._algorithm import AlgorithmParameters
from cal_disp.product import CalProduct, DispProduct

logger = logging.getLogger(__name__)


# Internal helpers
def _build_event_mask(
    disp_file: Path,
    defo_area_db: Path | None,
    event_db: Path | None,
    mask_dir: Path,
) -> Path | None:
    """Build one mask GeoTIFF (1 = valid) from the GeoJSON databases.

    Returns ``None`` if no database is given.
    """
    from cal_disp.prep.generate_event_mask import generate_event_mask

    mask_dir.mkdir(parents=True, exist_ok=True)
    mask_paths: list[Path] = []

    for db_path, label, filter_by_date in [
        (defo_area_db, "defo_area", False),  # always mask for the frame
        (event_db, "event", True),  # mask only within the epoch
    ]:
        if db_path is not None:
            out = generate_event_mask(
                product_path=disp_file,
                geojson_path=db_path,
                output_path=mask_dir / f"{label}_mask.tif",
                filter_by_date=filter_by_date,
            )
            mask_paths.append(out)

    if not mask_paths:
        return None
    if len(mask_paths) == 1:
        return mask_paths[0]

    combined = mask_dir / "combined_event_mask.tif"
    with rasterio.open(mask_paths[0]) as src:
        combined_mask = src.read(1).astype(bool)
        meta = src.meta.copy()
    for p in mask_paths[1:]:
        with rasterio.open(p) as src:
            combined_mask &= src.read(1).astype(bool)
    with rasterio.open(combined, "w", **meta) as dst:
        dst.write(combined_mask.astype(np.uint8)[np.newaxis])
    return combined


def _read_wavelength_m(disp_file: Path) -> float:
    """Read the radar wavelength (m) from a DISP product."""
    from cal_disp.product._disp import read_disp_metadata

    try:
        wl_m = read_disp_metadata(
            disp_file, {"radar_wavelength": "/identification/radar_wavelength"}
        )["radar_wavelength"]
    except OSError as e:
        raise RuntimeError(f"Could not open {disp_file.name}") from e
    if wl_m is None:
        msg = (
            f"Could not read /identification/radar_wavelength from {disp_file.name}. "
            "Ensure the file is a valid DISP product with an /identification group."
        )
        raise RuntimeError(msg)
    logger.debug("Read radar_wavelength %.6f m from %s", wl_m, disp_file.name)
    return float(wl_m)


def _unwrap_cycle_length_m(wavelength_m: float) -> float:
    """LOS displacement of one unwrapping cycle (2π of phase), in metres.

    Repeat-pass InSAR measures the two-way path, so one phase cycle is
    ``wavelength / 2`` of LOS displacement (27.7 mm for Sentinel-1 C-band,
    not the 55.5 mm wavelength).  Venti's ``UnwrapCorrector`` rounds region
    offsets to integer multiples of the value it is given as ``wavelength``,
    so it must receive this cycle length, not the wavelength itself.
    """
    return wavelength_m / 2.0


def _load_displacement(ds_disp: xr.Dataset, block_rows: int = 512) -> np.ndarray:
    """Load the ``displacement`` layer once, as float32.

    The layer is float64 on disk (585 MB for a DISP-S1 frame).  Reading it
    with ``.values`` loads the whole float64 array and keeps it cached in the
    dataset for the rest of the run, next to the float32 copy the fit uses.
    Row blocks are read into a preallocated float32 array instead, so the
    full float64 array is never held and nothing stays cached.  The values
    equal ``displacement.values.astype(np.float32)`` exactly.
    """
    return _load_layer(ds_disp, "displacement", np.float32, block_rows)


def _load_layer(
    ds_disp: xr.Dataset, name: str, dtype: type, block_rows: int = 512
) -> np.ndarray:
    """Load one 2-D layer in row blocks into a new `dtype` array.

    ``ds[name].values`` keeps the decoded array (often float32 or float64,
    even for masks) cached in the dataset for the rest of the run; reading
    row blocks into a preallocated array holds only the output. The values
    equal ``ds[name].values.astype(dtype)``.
    """
    layer = ds_disp[name]
    out = np.empty(layer.shape, dtype=dtype)
    for row in range(0, layer.shape[0], block_rows):
        out[row : row + block_rows] = layer[row : row + block_rows].values
    return out


def _remove_restore_valid(
    defo_area_db_json: Path | None,
    event_db_json: Path | None,
    x: np.ndarray,
    y: np.ndarray,
    utm_epsg: int,
    reference_date: Any,
    secondary_date: Any,
) -> np.ndarray | None:
    """Pixels outside the curated areas and the pair's events (True = valid).

    The v0.5 databases (Venti T32/T41: `id`/`name`/`version`, events with
    `t0`/`t1`) are read with Venti's loaders and its `Event.affects` rule,
    and rasterized on the product grid: cheaper than a point-in-polygon test
    on every pixel centre. None when neither database is given.
    """
    if defo_area_db_json is None and event_db_json is None:
        return None
    import shapely
    from pyproj import Transformer
    from rasterio.features import rasterize
    from rasterio.transform import from_origin
    from venti.calibration.remove_restore import load_area_db, load_event_db

    shapes = []
    if defo_area_db_json is not None:
        areas = load_area_db(defo_area_db_json)
        shapes += areas.shapes()
        logger.info("remove-restore areas: %s (%d)", areas.version, len(areas.items))
    if event_db_json is not None:
        events = load_event_db(event_db_json)
        active = events.active(reference_date, secondary_date)
        shapes += [e.shape() for e in active]
        logger.info(
            "remove-restore events: %s, %d of %d affect the pair %s",
            events.version,
            len(active),
            len(events.items),
            [e.id for e in active],
        )
    shape = (len(y), len(x))
    if not shapes:
        return np.ones(shape, dtype=bool)
    to_utm = Transformer.from_crs("EPSG:4326", f"EPSG:{utm_epsg}", always_xy=True)

    def _to_utm(lonlat: np.ndarray) -> np.ndarray:
        # shapely.transform hands over an (N, 2) coordinate array
        return np.column_stack(to_utm.transform(lonlat[:, 0], lonlat[:, 1]))

    dx, dy = float(x[1] - x[0]), float(y[0] - y[1])
    transform = from_origin(x[0] - dx / 2, y[0] + dy / 2, dx, dy)
    inside = rasterize(
        [(shapely.transform(g, _to_utm), 1) for g in shapes],
        out_shape=shape,
        transform=transform,
        fill=0,
        dtype="uint8",
    ).astype(bool)
    return ~inside


def _read_external_mask(mask_file: Path, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Read the runconfig ``mask_file``: True where the pixel may be used.

    The file must be a single-band uint8 raster on the DISP grid (same shape,
    pixel size and origin; 0 = invalid / no data, 1 = good). Anything else is
    an input error rather than a silent resampling, so the PGE notices a mask
    produced for another grid.
    """
    with rasterio.open(mask_file) as src:
        if src.count != 1:
            raise ValueError(
                f"mask_file {mask_file.name}: expected 1 band, got {src.count}"
            )
        if src.dtypes[0] != "uint8":
            raise ValueError(
                f"mask_file {mask_file.name}: dtype must be uint8, got {src.dtypes[0]}"
            )
        if src.shape != (len(y), len(x)):
            raise ValueError(
                f"mask_file {mask_file.name}: shape {src.shape} does not match the"
                f" DISP grid {(len(y), len(x))}"
            )
        dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
        t = src.transform
        # pixel-corner origin of the DISP grid (x/y are pixel centres)
        x0, y0 = float(x[0]) - dx / 2, float(y[0]) - dy / 2
        tol = 0.01 * abs(dx)
        if (
            abs(t.a - dx) > tol
            or abs(t.e - dy) > tol
            or abs(t.c - x0) > tol
            or abs(t.f - y0) > tol
        ):
            raise ValueError(
                f"mask_file {mask_file.name}: geotransform {tuple(t)[:6]} is not"
                f" the DISP grid (origin {x0}, {y0}; spacing {dx}, {dy})"
            )
        data = src.read(1)
    bad = ~np.isin(data, (0, 1))
    if bad.any():
        raise ValueError(
            f"mask_file {mask_file.name}: values must be 0 or 1, found"
            f" {np.unique(data[bad])[:5].tolist()}"
        )
    return data == 1


def _read_dem_on_grid(dem_file: Path, shape: tuple[int, int]) -> np.ndarray:
    """Read the DISP-S1-STATIC DEM (float64, NaN for nodata); must be on the grid."""
    with rasterio.open(dem_file) as src:
        if src.shape != tuple(shape):
            raise ValueError(
                f"dem_file {Path(dem_file).name}: shape {src.shape} does not match"
                f" the DISP grid {tuple(shape)}"
            )
        dem = src.read(1, masked=True).astype(np.float64)
    return np.ma.filled(dem, np.nan)


def _date_to_decimal_year(dt: datetime) -> float:
    """Convert a datetime to a decimal year (e.g. 2022.55)."""
    year = dt.year
    year_start = datetime(year, 1, 1, tzinfo=dt.tzinfo)
    year_end = datetime(year + 1, 1, 1, tzinfo=dt.tzinfo)
    elapsed = (dt - year_start).total_seconds()
    total = (year_end - year_start).total_seconds()
    return year + elapsed / total


def _load_los_bands(
    los_file: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read the (east, north, up) bands of a LOS GeoTIFF; NaN where no data."""
    with rasterio.open(los_file) as src:
        if src.count < 3:
            raise ValueError(
                "LOS file must have ≥3 bands (east, north, up),"
                f" got {src.count}: {los_file}"
            )
        los_east = src.read(1).astype(np.float32)
        los_north = src.read(2).astype(np.float32)
        los_up = src.read(3).astype(np.float32)
        nodata = src.nodata
    from cal_disp.prep.consistency import mask_los_nodata

    return mask_los_nodata(los_east, los_north, los_up, nodata)


def _staged_station_files(
    unr_timeseries_dir: Path, reference_frame: str, grid_type: str
) -> list[Path]:
    """Staged UNR files of one reference frame and grid type."""
    return sorted(unr_timeseries_dir.glob(f"*_{reference_frame}_{grid_type}.tenv8"))


def _gnss_los_fields(
    gnss_ref,
    los: tuple[np.ndarray, np.ndarray, np.ndarray],
    disp_file: Path,
    x: np.ndarray,
    y: np.ndarray,
    factor: int,
    cache_dir: Path,
    ref_date: float,
    sec_date: float,
) -> tuple[np.ndarray, np.ndarray]:
    """GNSS LOS displacement and its uncertainty (m) on the product grid.

    With ``factor > 1`` both are computed at the block centres of the fit grid
    and interpolated back: the fit only sees block averages anyway, and it is
    much faster (NYC F08622: <= 0.01 mm change, ~93 min -> seconds).
    """
    from venti.gnss import compute_gnss_los, compute_gnss_los_std

    los_east, los_north, los_up = los
    grid: Path | tuple[np.ndarray, np.ndarray] = disp_file
    c = factor // 2
    if factor > 1:
        grid = (x[c::factor], y[c::factor])
        los_east, los_north, los_up = (a[c::factor, c::factor] for a in los)
    kwargs = {
        "gnss_ref": gnss_ref,
        "los_east": los_east,
        "los_north": los_north,
        "los_up": los_up,
        "grid": grid,
        "cache_dir": cache_dir,
        "ref_date": ref_date,
        "sec_date": sec_date,
    }
    fields = [compute_gnss_los(**kwargs), compute_gnss_los_std(**kwargs)]
    if factor > 1:
        fields = [
            _upsample_centres(f, x[c::factor], y[c::factor], x, y) for f in fields
        ]
    gnss_los, gnss_los_std = (np.asarray(f, dtype=np.float32) / 1000.0 for f in fields)
    return gnss_los, gnss_los_std


def _gnss_los_fields_enu(
    gnss_options: Any,
    grid_lookup: Path,
    station_dir: Path,
    utm_epsg: int,
    reference_frame: str,
    grid_version: str,
    los: tuple[np.ndarray, np.ndarray, np.ndarray],
    x: np.ndarray,
    y: np.ndarray,
    factor: int,
    dt_years: float,
    defo_area_db_json: Path | None = None,
) -> tuple[np.ndarray, np.ndarray, Any]:
    """GNSS LOS displacement and sigma (m) through Venti's single sampling path.

    `venti.gnss.sampling.sample_gnss_enu` on the fit-grid block centres: grid
    nodes within ``gnss.buffer_meters`` outside the frame constrain the
    edges (R-G2), E/N/U are interpolated first and projected per pixel, so
    nodes outside the swath need no LOS look, and nodes inside the curated
    defo areas can be dropped or re-interpolated (R-G5). Returns the two
    fields and the provenance record.
    """
    from venti.gnss.sampling import (
        GnssGridConfig,
        project_field_to_los,
        sample_gnss_enu,
    )

    c = factor // 2
    los_east, los_north, los_up = los
    gx, gy = x, y
    if factor > 1:
        gx, gy = x[c::factor], y[c::factor]
        los_east, los_north, los_up = (a[c::factor, c::factor] for a in los)
    cfg = GnssGridConfig(
        grid_lookup=grid_lookup,
        station_dir=station_dir,
        utm_epsg=utm_epsg,
        reference_frame=reference_frame,
        grid_type="constant",
        version=grid_version,
        buffer_meters=gnss_options.buffer_meters,
        exclude_defo_nodes=gnss_options.exclude_defo_nodes,
        reinterpolate_excluded=gnss_options.reinterpolate_excluded,
        snapshot_id=gnss_options.snapshot_id,
    )
    exclude = db_version = None
    if gnss_options.exclude_defo_nodes and defo_area_db_json is not None:
        from venti.calibration.remove_restore import load_area_db

        db = load_area_db(defo_area_db_json)
        exclude, db_version = db.shapes(), db.version
    field = sample_gnss_enu(cfg, (gx, gy), exclude, db_version)
    los_mm, sigma_mm = project_field_to_los(
        field, los_east, los_north, los_up, dt_years=dt_years
    )
    fields = [los_mm, sigma_mm]
    if factor > 1:
        fields = [_upsample_centres(f, gx, gy, x, y) for f in fields]
    gnss_los, gnss_los_std = (np.asarray(f, dtype=np.float32) / 1000.0 for f in fields)
    return gnss_los, gnss_los_std, field.provenance


#: Pixel step of the coarse `/auxiliary` grid (about 5 km at 30 m posting).
AUXILIARY_STEP = 167


def _auxiliary_gnss_displacement(
    gnss_options: Any,
    grid_lookup: Path,
    station_dir: Path,
    utm_epsg: int,
    reference_frame: str,
    grid_version: str,
    x: np.ndarray,
    y: np.ndarray,
    dt_years: float,
    defo_area_db_json: Path | None = None,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """GNSS E/N/U displacement over the pair and its sigma (m) on the coarse grid.

    The UNR constant grid sampled by `venti.gnss.sampling.sample_gnss_enu` with
    the calibration's own options (buffer, defo-node exclusion and
    re-interpolation), at the `/auxiliary` grid points, and scaled by the pair
    span exactly as `project_field_to_los` scales the LOS field: velocity x
    ``dt_years`` and sigma x ``|dt_years|``. Divide by the span for velocity.
    """
    from venti.gnss.sampling import GnssGridConfig, sample_gnss_enu

    cfg = GnssGridConfig(
        grid_lookup=grid_lookup,
        station_dir=station_dir,
        utm_epsg=utm_epsg,
        reference_frame=reference_frame,
        grid_type="constant",
        version=grid_version,
        buffer_meters=gnss_options.buffer_meters,
        exclude_defo_nodes=gnss_options.exclude_defo_nodes,
        reinterpolate_excluded=gnss_options.reinterpolate_excluded,
        snapshot_id=gnss_options.snapshot_id,
    )
    exclude = db_version = None
    if gnss_options.exclude_defo_nodes and defo_area_db_json is not None:
        from venti.calibration.remove_restore import load_area_db

        db = load_area_db(defo_area_db_json)
        exclude, db_version = db.shapes(), db.version
    field = sample_gnss_enu(
        cfg, (x[::AUXILIARY_STEP], y[::AUXILIARY_STEP]), exclude, db_version
    )
    names = {"east_west": "e", "north_south": "n", "up_down": "u"}
    disp = {
        k: (np.asarray(getattr(field, f"v{c}")) * dt_years / 1000.0).astype(np.float32)
        for k, c in names.items()
    }
    std = {
        f"{k}_std": (
            (np.asarray(getattr(field, f"sigma_v{c}")) * abs(dt_years) / 1000.0).astype(
                np.float32
            )
        )
        for k, c in names.items()
    }
    return disp, std


def _plate_motion_los(
    plate: str,
    utm_epsg: int,
    los: Path | tuple[np.ndarray, np.ndarray, np.ndarray],
    x: np.ndarray,
    y: np.ndarray,
    dt_years: float,
    block_rows: int = 512,
) -> np.ndarray:
    """LOS displacement (m) of rigid `plate` motion (ITRF2020-PMM) over the pair.

    The plate's ENU velocity (`geepers.euler.plate_velocity_enu`, mm/yr, no up
    component) is evaluated on a ~5 km grid spanning the frame, where a rigid
    rotation is linear to well below a micron, interpolated to each pixel and
    projected on its LOS, row block by row block to keep memory flat.
    ``DISP - calibration - plate_motion`` is the plate-fixed displacement.
    NaN where the LOS has no look.

    `los` is the LOS GeoTIFF (read one row window at a time, so no full LOS
    frame is held: this runs after the fit, under the peak-memory budget) or
    the three bands already in memory.
    """
    from geepers.euler import plate_velocity_enu
    from pyproj import Transformer

    # uniform coarse grid spanning the frame edge to edge, so no pixel is
    # beyond the outermost sample (interpolation, never edge-hold)
    xc = np.linspace(x[0], x[-1], -(-len(x) // AUXILIARY_STEP) + 1)
    yc = np.linspace(y[0], y[-1], -(-len(y) // AUXILIARY_STEP) + 1)
    xx, yy = np.meshgrid(xc, yc)
    lon, lat = Transformer.from_crs(utm_epsg, 4326, always_xy=True).transform(xx, yy)
    ve, vn, _ = plate_velocity_enu(lon.ravel(), lat.ravel(), plate)
    ve, vn = (np.asarray(v, dtype=np.float64).reshape(xx.shape) for v in (ve, vn))
    out = np.empty((len(y), len(x)), dtype=np.float32)
    scale = dt_years / 1000.0
    bands = None if isinstance(los, (str, Path)) else los
    src = rasterio.open(los) if bands is None else None
    try:
        for r0 in range(0, len(y), block_rows):
            rows = slice(r0, min(r0 + block_rows, len(y)))
            if src is not None:
                window = Window(0, rows.start, len(x), rows.stop - rows.start)
                e, n, u = src.read((1, 2, 3), window=window).astype(np.float32)
                no_look = (e == 0) & (n == 0) & (u == 0)
                if src.nodata is not None and not np.isnan(src.nodata):
                    no_look |= (e == src.nodata) | (n == src.nodata) | (u == src.nodata)
            else:
                assert bands is not None
                e, n = bands[0][rows], bands[1][rows]
                no_look = (e == 0) & (n == 0)
            ve_b = _upsample_centres(ve, xc, yc, x, y[rows])
            vn_b = _upsample_centres(vn, xc, yc, x, y[rows])
            block = (e * ve_b + n * vn_b) * scale
            block[no_look] = np.nan
            out[rows] = block
    finally:
        if src is not None:
            src.close()
    # 1 micron: the smooth field's float32 mantissa noise otherwise defeats the
    # compression (golden pair: 35 MB -> 6 MB, error <= 0.5 micron)
    np.round(out, 6, out=out)
    return out


def _upsample_centres(
    field: np.ndarray, xc: np.ndarray, yc: np.ndarray, x: np.ndarray, y: np.ndarray
) -> np.ndarray:
    """Bilinear interpolation from grid (`yc`, `xc`) to (`y`, `x`).

    Beyond the outermost samples the edge value is kept.
    """

    def _weights(src: np.ndarray, dst: np.ndarray):
        n = len(src)
        if n == 1:
            zeros = np.zeros(len(dst), dtype=int)
            return zeros, zeros, np.zeros(len(dst), dtype=np.float32)
        pos = np.clip((dst - src[0]) / (src[1] - src[0]), 0, n - 1)
        i0 = np.minimum(pos.astype(int), n - 2)
        return i0, i0 + 1, (pos - i0).astype(np.float32)

    field = np.asarray(field, dtype=np.float32)
    x0, x1, wx = _weights(xc, x)
    y0, y1, wy = _weights(yc, y)
    rows = field[:, x0] * (1 - wx) + field[:, x1] * wx
    return rows[y0] * (1 - wy)[:, None] + rows[y1] * wy[:, None]


def _quality_reference_point(
    quality: np.ndarray,
    max_quality_test: float = 0.95,
    min_quality: float = 0.5,
    step: float = 0.05,
) -> tuple[int, int]:
    """Pick a high-quality pixel near the centre of mass of the good data.

    Port of ``opera_utils.disp.rebase_reference.find_reference_point`` (same
    thresholds and rule), working on the array instead of a GeoTIFF so the
    calibration core does not need ``opera-utils[disp]`` (dask, zarr, s3fs).
    Thresholds step down from `max_quality_test`; at the first one with
    candidates, the candidate closest to their centre of mass wins.

    Raises
    ------
    ValueError
        If no pixel exceeds `min_quality`.

    """
    from scipy import ndimage

    for threshold in np.arange(max_quality_test, min_quality - step, -step):
        candidates = quality > threshold
        if not np.any(candidates):
            continue
        r, c = ndimage.center_of_mass(candidates)
        rows, cols = np.where(candidates)
        idx = np.argmin(np.sqrt((rows - r) ** 2 + (cols - c) ** 2))
        return int(rows[idx]), int(cols[idx])
    msg = f"No valid candidates found with quality above {min_quality}"
    raise ValueError(msg)


def _build_unwrap_hook(
    ds_disp: xr.Dataset,
    disp_2d: np.ndarray,
    gnss_los: np.ndarray,
    mask: np.ndarray,
    venti_options: Any,
    posting_meters: float,
    *,
    pair: str = "",
) -> Any:
    """Build the whole-cycle unwrap hook of the loclin route (Venti T36, R-U1).

    Regions come from the DISP ``water_mask`` (land without
    ``temporal_coherence`` does not join regions); the medians use the same
    valid mask as the fit. ``timeseries_inversion_residuals`` is read only when
    the residual gate is configured.
    """
    from venti.unwrap import unwrap_hook_for_pair

    land = (
        ds_disp.water_mask.values.astype(bool)
        if "water_mask" in ds_disp
        else np.ones(disp_2d.shape, dtype=bool)
    )
    data_valid = (
        np.isfinite(ds_disp.temporal_coherence.values)
        if "temporal_coherence" in ds_disp
        else np.isfinite(disp_2d)
    )
    residual = None
    if venti_options.unwrap.residual_gate_cycles is not None:
        if "timeseries_inversion_residuals" not in ds_disp:
            msg = (
                "unwrap.residual_gate_cycles is set but the DISP product has no"
                " timeseries_inversion_residuals"
            )
            raise ValueError(msg)
        residual = ds_disp.timeseries_inversion_residuals.values
    return unwrap_hook_for_pair(
        np.squeeze(land),
        np.squeeze(data_valid),
        gnss_los,
        mask,
        venti_options,
        posting_meters,
        inversion_residual=np.squeeze(residual) if residual is not None else None,
        pair=pair,
    )


def _find_reference_point(
    coherence: np.ndarray | None, mask: np.ndarray
) -> tuple[int, int]:
    """Choose the reference pixel ``(row, col)`` for one product.

    Venti's rule (`_quality_reference_point`, ported from ``opera_utils``)
    applied to the product's temporal coherence over valid pixels (all
    pixels equal when the product has no coherence layer).
    """
    if coherence is not None:
        quality = np.asarray(coherence, dtype=np.float32)
    else:
        quality = np.ones(mask.shape, dtype=np.float32)
    quality = np.where(mask & np.isfinite(quality), quality, 0.0).astype(np.float32)
    try:
        row, col = _quality_reference_point(quality)
    except ValueError:
        # Nothing above the minimum quality: take the best pixel
        row, col = np.unravel_index(np.argmax(quality), quality.shape)
    ref_point = (int(row), int(col))
    logger.info("Reference pixel (row, col): %s", ref_point)
    return ref_point


def _setup_gnss_reference(
    disp_product: DispProduct,
    unr_grid_latlon_file: Path,
    unr_timeseries_dir: Path,
    gnss_dir: Path,
    reference_frame: str,
    grid_type: str,
):
    """Create a ``GNSSReference`` that uses the pre-staged UNR files.

    The staged files are symlinked into `gnss_dir`, where Venti's
    ``download_stations()`` finds them instead of downloading.
    """
    from venti.gnss import GNSSReference

    gnss_dir.mkdir(parents=True, exist_ok=True)

    # Venti's expected lookup name
    lookup_link = gnss_dir / "grid_latlon_lookup.txt"
    if not lookup_link.exists() and unr_grid_latlon_file.exists():
        lookup_link.symlink_to(unr_grid_latlon_file.resolve())

    staged = _staged_station_files(unr_timeseries_dir, reference_frame, grid_type)
    if not staged:
        msg = (
            f"No staged UNR files '*_{reference_frame}_{grid_type}.tenv8' in"
            f" {unr_timeseries_dir}. Stage them with `cal-disp download unr"
            f" --grid-type {grid_type}`."
        )
        raise FileNotFoundError(msg)
    for tenv8 in staged:
        link = gnss_dir / tenv8.name
        if not link.exists():
            link.symlink_to(tenv8.resolve())

    # UTM bounds (south, north, west, east)
    b = disp_product.get_bounds()
    bounds = (b["bottom"], b["top"], b["left"], b["right"])
    utm_epsg = disp_product.get_epsg()

    gnss_ref = GNSSReference(
        bounds=bounds,
        output_dir=gnss_dir,
        reference_frame=reference_frame,
        utm_epsg=utm_epsg,
        grid_type=grid_type,
    )
    n_stations = gnss_ref.download_stations()
    logger.info("GNSS setup complete: %d stations available", n_stations)
    return gnss_ref


def _in_scratch_temp_dir(func):
    """Run `func` with all temp files (Python and GDAL) in ``<work_directory>/tmp``.

    The system temp dir or cwd may not be writable in a PGE container.
    """
    signature = inspect.signature(func)

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        from venti.workflow.utils import scratch_temp_dir

        bound = signature.bind(*args, **kwargs)
        bound.apply_defaults()
        work_directory = bound.arguments["work_directory"]
        if work_directory is None:
            work_directory = Path(bound.arguments["output_dir"]) / "scratch"
        with scratch_temp_dir(work_directory):
            return func(*args, **kwargs)

    return wrapper


# Public entry point
@_in_scratch_temp_dir
def run_calibration(
    disp_file: Path,
    unr_grid_latlon_file: Path,
    unr_timeseries_dir: Path,
    output_dir: Path,
    algorithm_parameters: AlgorithmParameters | None = None,
    dem_file: Path | None = None,
    los_file: Path | None = None,
    reference_tropo_files: list[Path] | None = None,
    secondary_tropo_files: list[Path] | None = None,
    defo_area_db_json: Path | None = None,
    event_db_json: Path | None = None,
    fit_mask_file: Path | None = None,
    mask_file: Path | None = None,
    block_shape: tuple[int, int] = (512, 512),  # noqa: ARG001 — TODO: Dask blocking
    n_workers: int = 4,  # noqa: ARG001 — TODO: Dask parallelisation
    threads_per_worker: int = 1,
    work_directory: Path | None = None,
    pge_runconfig: str | None = None,
    frame_parameters_version: str | None = None,
    # Calibration reference metadata
    calibration_reference_name: str = "UNR gridded data",
    calibration_reference_version: str = "0.3",
    calibration_reference_type: str | None = None,
    calibration_reference_reference_frame: str = "IGS20",
    # Product metadata. Platform, orbit, track, look direction, instrument,
    # band, DEM, imaging geometry, satellite names and CEOS fields are read
    # from the DISP product's identification group.
    product_version: str = "1.0",
    compression: bool = True,
    processing_facility: str = "NASA Jet Propulsion Laboratory on AWS",
    product_data_access: str = "https://search.asf.alaska.edu/#/?dataset=OPERA-S1&productTypes=DISP-S1-CAL",
    static_layers_data_access: str | None = None,
    source_data_access: str | None = None,
) -> Path:
    """Run the single-file displacement calibration workflow.

    Parameters
    ----------
    disp_file : Path
        Input OPERA L3 DISP-S1 product (NetCDF).
    unr_grid_latlon_file : Path
        UNR grid station lookup table (``grid_latlon_lookup.txt``).
    unr_timeseries_dir : Path
        Directory containing pre-staged UNR ``.tenv8`` station files.
    output_dir : Path
        Directory for the output ``CalProduct`` NetCDF.
    algorithm_parameters : AlgorithmParameters, optional
        Algorithm configuration.  Defaults are used when ``None``.
    dem_file : Path, optional
        DEM GeoTIFF; required with tropo files.
    los_file : Path
        LOS GeoTIFF (bands: east, north, up).
    reference_tropo_files, secondary_tropo_files : list[Path], optional
        OPERA TROPO-ZENITH files (1-2 per date). With both, the differential
        tropospheric delay is removed before the fit (unless
        ``apply_tropo_correction`` is false).
    defo_area_db_json : Path, optional
        GeoJSON of deformation areas to exclude from the fit.
    event_db_json : Path, optional
        GeoJSON of events to exclude from the fit (within the epoch only).
    fit_mask_file : Path, optional
        Stable-pixel fit mask (uint8 GeoTIFF on the DISP grid, 1 = may be fit)
        from the runconfig, for masks not packaged in Venti; loclin only and
        exclusive with ``surface.fit_mask``.
    mask_file : Path, optional
        Runconfig byte mask (uint8, 0 = invalid, 1 = good) on the DISP grid;
        its 0 pixels are excluded from the fit like the product's own masks.
    block_shape : tuple[int, int]
        Unused (reserved for Dask).
    n_workers : int
        Unused (reserved for Dask).
    threads_per_worker : int
        Parallel jobs for the windowed fit.
    work_directory : Path, optional
        Scratch directory, by default ``output_dir / "scratch"``.
    pge_runconfig : str, optional
        Serialised PGE run-config YAML stored in product metadata.
    frame_parameters_version : str, optional
        Version of the frame-parameter table applied (the overrides file's
        top-level ``version``); recorded in ``/metadata`` on the v0.5 route.
    calibration_reference_name : str
        Human-readable name of the calibration reference dataset.
    calibration_reference_version : str
        Version string of the calibration reference dataset.
    calibration_reference_type : str, optional
        ``'constant'`` or ``'variable'``; must match ``grid_type`` (default).
    calibration_reference_reference_frame : str
        GNSS reference frame (e.g. ``'IGS20'``).
    product_version : str
        Product version in ``<major>.<minor>`` format; written to the product
        filename and identification metadata.
    compression : bool
        Write the product rasters gzip-compressed in (256, 256) chunks.
    processing_facility : str
        Processing facility name.
    product_data_access : str
        URL for product data access.
    static_layers_data_access : str, optional
        URL of the frame's static layers; by default the DISP product's value.
    source_data_access : str, optional
        URL for source (DISP) data access; by default the DISP product's
        ``product_data_access``.

    Returns
    -------
    Path
        Path to the output ``CalProduct`` NetCDF file.

    """
    if los_file is None:
        raise ValueError("los_file is required for calibration")

    if algorithm_parameters is None:
        algorithm_parameters = AlgorithmParameters()

    if work_directory is None:
        work_directory = output_dir / "scratch"
    work_directory.mkdir(parents=True, exist_ok=True)
    processing_start = datetime.now(tz=timezone.utc)

    cal = algorithm_parameters.calibration_options
    if calibration_reference_type is None:
        calibration_reference_type = cal.grid_type
    elif calibration_reference_type != cal.grid_type:
        msg = (
            f"UNR data type '{calibration_reference_type}' does not match the"
            f" calibration grid_type '{cal.grid_type}'"
        )
        raise ValueError(msg)

    # Load DISP product
    logger.info("Loading DISP product: %s", disp_file.name)
    disp_product = DispProduct.from_path(disp_file)
    # The static layers are used by array index: same grid and frame required
    from cal_disp.prep.consistency import check_input_geometry

    check_input_geometry(disp_product, los_file, dem_file)
    ds_disp = disp_product.open_dataset()

    time = ds_disp.time.values
    y = ds_disp.y.values
    x = ds_disp.x.values

    spatial_ref = ds_disp.get("spatial_ref")

    x_spacing = float(np.abs(x[1] - x[0]))
    y_spacing = float(np.abs(y[1] - y[0]))
    product_sample_spacing = f"{x_spacing:g}m"

    # Grid extent: outer pixel edges in UTM, and the same in lon/lat (as DISP)
    from cal_disp.product.output._utils import bounding_polygon_wkt, grid_bounds

    west, south, east, north = grid_bounds(x, y)
    product_bounding_box = f"({west}, {south}, {east}, {north})"
    crs_wkt = spatial_ref.attrs.get("crs_wkt") if spatial_ref is not None else None
    if crs_wkt:
        bounding_polygon = bounding_polygon_wkt(x, y, crs_wkt)
    else:
        logger.warning("%s has no CRS; bounding_polygon left empty", disp_file.name)
        bounding_polygon = ""

    source_data_file_list = [disp_file.name]
    source_calibration_file_list = [unr_grid_latlon_file.name]
    tenv8_files = _staged_station_files(
        unr_timeseries_dir, cal.reference_frame, cal.grid_type
    )
    source_calibration_file_list.extend(f.name for f in tenv8_files)

    # Platform/orbit/geometry metadata from the DISP product (WARNING + marked
    # fallback for anything the input lacks)
    disp_meta = disp_product.read_metadata()

    def _meta(name: str, fallback: str | int) -> str | int:
        value = disp_meta.get(name)
        return fallback if value is None else value

    platform_id = str(_meta("platform_id", "unknown"))
    satellite_names = str(_meta("source_data_satellite_names", platform_id))
    source_data_satellite_names = [s.strip() for s in satellite_names.split(",")]
    if static_layers_data_access is None:
        static_layers_data_access = str(_meta("static_layers_data_access", "unknown"))
    if source_data_access is None:
        source_data_access = str(_meta("product_data_access", "unknown"))

    # GNSS reference setup
    logger.info("Setting up GNSS reference...")
    gnss_dir = work_directory / "gnss"
    gnss_ref = _setup_gnss_reference(
        disp_product=disp_product,
        unr_grid_latlon_file=unr_grid_latlon_file,
        unr_timeseries_dir=unr_timeseries_dir,
        gnss_dir=gnss_dir,
        reference_frame=cal.reference_frame,
        grid_type=cal.grid_type,
    )

    # Load LOS unit vectors
    logger.info("Loading LOS unit vectors: %s", los_file.name)
    los_east, los_north, los_up = _load_los_bands(los_file)

    # GNSS LOS reference for this interval (Venti returns mm; disp is in m)
    from venti.calibration.tropo import apply_tropo, choose_tropo_mode
    from venti.calibration.two_pass import calibrate_pair
    from venti.io import read_netcdf_correction
    from venti.surface import SENTINEL1_WAVELENGTH_M

    ref_decimal = _date_to_decimal_year(disp_product.reference_date)
    sec_decimal = _date_to_decimal_year(disp_product.secondary_date)
    v_opts = cal.to_venti()
    gnss_provenance = None  # the gamma route has no provenance record
    if v_opts.surface.method == "loclin" and cal.grid_type == "constant":
        # v0.5: the single sampling path with the grid buffer (R-G2, R-G5)
        gnss_los, gnss_los_std, gnss_provenance = _gnss_los_fields_enu(
            v_opts.gnss,
            Path(unr_grid_latlon_file),
            gnss_dir,
            disp_product.get_epsg(),
            cal.reference_frame,
            calibration_reference_version,
            (los_east, los_north, los_up),
            x,
            y,
            cal.downsample_factor,
            sec_decimal - ref_decimal,
            defo_area_db_json=defo_area_db_json,
        )
        logger.info(
            "GNSS field: %d grid nodes, buffer %.0f m, %d excluded / %d"
            " re-interpolated (provenance %s)",
            gnss_provenance.n_nodes,
            gnss_provenance.buffer_meters,
            gnss_provenance.n_excluded,
            gnss_provenance.n_reinterpolated,
            gnss_provenance.digest[:12],
        )
    else:
        # gamma route: unchanged, so the gamma golden stays bit-identical
        gnss_los, gnss_los_std = _gnss_los_fields(
            gnss_ref,
            (los_east, los_north, los_up),
            disp_file,
            x,
            y,
            cal.downsample_factor,
            gnss_dir,
            ref_decimal,
            sec_decimal,
        )
    # /auxiliary: the GNSS E/N/U displacement over the pair, on the coarse grid
    # (constant grid only: the variable grid is not sampled through Venti)
    model_3d = model_3d_std = None
    if cal.grid_type == "constant":
        model_3d, model_3d_std = _auxiliary_gnss_displacement(
            v_opts.gnss,
            Path(unr_grid_latlon_file),
            gnss_dir,
            disp_product.get_epsg(),
            cal.reference_frame,
            calibration_reference_version,
            x,
            y,
            sec_decimal - ref_decimal,
            defo_area_db_json=defo_area_db_json,
        )
    # the LOS bands only serve the GNSS projection; three full float32 frames
    del los_east, los_north, los_up

    # Displacement (m): float64 on disk, read in row blocks straight into
    # float32 so the float64 array is never held or cached in `ds_disp`
    disp_2d = _load_displacement(ds_disp)

    # Valid pixels (both masks: 1 = valid)
    mask = ~np.isnan(disp_2d)
    if "recommended_mask" in ds_disp:
        mask &= _load_layer(ds_disp, "recommended_mask", bool)
    if "water_mask" in ds_disp:
        mask &= _load_layer(ds_disp, "water_mask", bool)
    if mask_file is not None:
        external = _read_external_mask(Path(mask_file), x, y)
        logger.info(
            "mask_file %s: %d pixels marked invalid",
            Path(mask_file).name,
            int((~external).sum()),
        )
        mask &= external

    # Event areas (True = valid) are filled from neighbours before the fit
    event_mask = None
    if cal.to_venti().surface.method == "loclin":
        # v0.5: the Venti databases (T32/T41) and Venti's event rule
        event_mask = _remove_restore_valid(
            defo_area_db_json,
            event_db_json,
            x,
            y,
            disp_product.get_epsg(),
            disp_product.reference_date,
            disp_product.secondary_date,
        )
        if event_mask is not None:
            logger.info("Remove-restore: %d pixels to fill", int((~event_mask).sum()))
    elif defo_area_db_json is not None or event_db_json is not None:
        # gamma route: the gamma-era GeoJSON (frame_id / event_date), unchanged
        event_mask_file = _build_event_mask(
            disp_file=disp_file,
            defo_area_db=defo_area_db_json,
            event_db=event_db_json,
            mask_dir=work_directory / "event_masks",
        )
        if event_mask_file is not None:
            with rasterio.open(event_mask_file) as _src:
                event_mask = _src.read(1).astype(bool)
            logger.info("Event mask: %d pixels to fill", int((~event_mask).sum()))

    # Signals GNSS does not contain: removed before the fit and restored into
    # the calibration as components (Venti calibrate_pair; calibration ==
    # sum of components)
    venti_options = cal.to_venti()
    from venti.workflow.config import ALGORITHM_SCHEMA_VERSION

    have_tropo_files = bool(reference_tropo_files and secondary_tropo_files)
    dem_2d = None
    if venti_options.tropo.mode in ("auto", "stratified") and have_tropo_files:
        if dem_file is None:
            raise ValueError(
                f"tropo mode {venti_options.tropo.mode!r} needs dem_file (DEM relief"
                " and the stratified fit)"
            )
        dem_2d = _read_dem_on_grid(dem_file, disp_2d.shape)
    tropo_mode = choose_tropo_mode(
        venti_options.tropo,
        dem_2d,
        apply_tropo_correction=cal.apply_tropo_correction,
        valid=mask,
    )
    tropo_component = None
    if tropo_mode != "off" and not have_tropo_files:
        logger.info("tropo mode %s but no tropo files; calibrating without", tropo_mode)
    elif tropo_mode == "off" and have_tropo_files:
        logger.info("tropo mode off; ignoring the tropo files")
    elif tropo_mode != "off":
        if dem_file is None:
            raise ValueError(
                "dem_file is required when tropospheric correction files are provided"
            )
        from cal_disp.prep.tropo import prepare_troposphere_correction

        logger.info(
            "Preparing tropospheric correction (ref + sec), mode %s...", tropo_mode
        )
        tropo_dir = work_directory / "troposphere"
        ref_tropo_path, sec_tropo_path = prepare_troposphere_correction(
            disp_file=disp_file,
            dem_file=dem_file,
            los_file=los_file,
            reference_tropo_files=reference_tropo_files,
            secondary_tropo_files=secondary_tropo_files,
            output_dir=tropo_dir,
        )
        with rasterio.open(ref_tropo_path) as _src:
            ref_tropo = _src.read(1).astype(np.float32)
        with rasterio.open(sec_tropo_path) as _src:
            sec_tropo = _src.read(1).astype(np.float32)
        tropo_component = apply_tropo(
            sec_tropo - ref_tropo, tropo_mode, dem=dem_2d, valid=mask, x=x, y=y
        )

    set_component = None
    if cal.apply_solid_earth_tide_correction:
        set_corr = read_netcdf_correction(disp_file, "solid_earth_tide")
        if set_corr is None:
            logger.warning(
                "No /corrections/solid_earth_tide in %s; calibrating without SET"
                " correction",
                disp_file.name,
            )
        else:
            set_component = set_corr.astype(np.float32, copy=False)
            del set_corr

    coherence = None
    if "temporal_coherence" in ds_disp:
        coherence = _load_layer(ds_disp, "temporal_coherence", np.float32)
    elif cal.downsample_factor > 1 and cal.downsample_weighted:
        raise ValueError(
            "downsample_weighted=True but temporal_coherence is not in"
            f" {disp_file.name}"
        )

    ref_point = _find_reference_point(coherence, mask)
    # Radar wavelength from the product so other sensors (e.g. NISAR) are
    # handled too; one unwrapping cycle is lambda/2 of LOS displacement
    # (`_unwrap_cycle_length_m`).
    wavelength_m = (
        _read_wavelength_m(disp_file)
        if cal.unwrap_error_correction
        else SENTINEL1_WAVELENGTH_M
    )

    logger.info(
        "Fitting calibration surface (%s, window=%d px, downsample=%d, tropo=%s,"
        " SET=%s)...",
        venti_options.surface.method,
        cal.window_size_pixels,
        cal.downsample_factor,
        tropo_mode if tropo_component is not None else "off",
        set_component is not None,
    )
    unwrap_hook = None
    if cal.unwrap_error_correction and venti_options.surface.method == "loclin":
        unwrap_hook = _build_unwrap_hook(
            ds_disp,
            disp_2d,
            gnss_los,
            mask,
            venti_options,
            cal.posting_meters,
            pair=disp_file.name,
        )
    # a frame's stable-pixel mask, packaged (Venti surface.fit_mask, frame
    # table) or from the runconfig (fit_mask_file): unstable pixels leave the
    # fit but keep the surface
    fit_mask = None
    if fit_mask_file is not None:
        if venti_options.surface.fit_mask is not None:
            raise ValueError(
                f"fit_mask_file {fit_mask_file} and surface.fit_mask"
                f" {venti_options.surface.fit_mask!r} are both set; give one"
            )
        if venti_options.surface.method != "loclin":
            raise ValueError(
                "fit_mask_file needs surface.method 'loclin' (got"
                f" {venti_options.surface.method!r})"
            )
        fit_mask = _read_external_mask(Path(fit_mask_file), x, y)
        logger.info(
            "fit_mask_file %s: %d of %d pixels outside it are left out of the fit",
            Path(fit_mask_file).name,
            int((~fit_mask & mask).sum()),
            int(mask.sum()),
        )
    elif venti_options.surface.fit_mask is not None:
        from venti.frames import fit_mask_path

        fit_mask = _read_external_mask(
            fit_mask_path(venti_options.surface.fit_mask), x, y
        )
        logger.info(
            "fit_mask %s: %d of %d pixels outside it are left out of the fit",
            venti_options.surface.fit_mask,
            int((~fit_mask & mask).sum()),
            int(mask.sum()),
        )
    try:
        result = calibrate_pair(
            disp_2d,
            gnss_los,
            mask,
            ref_point,
            venti_options,
            cal.posting_meters,
            _unwrap_cycle_length_m(wavelength_m),
            gnss_los_std=gnss_los_std,
            coherence=coherence,
            tropo=tropo_component,
            set_correction=set_component,
            exclude_mask=~event_mask if event_mask is not None else None,
            fit_mask=fit_mask,
            n_jobs=threads_per_worker,
            unwrap_hook=unwrap_hook,
            # the product carries calibration and sigma_CAL only
            diagnostic_maps=False,
        )
    except ValueError as exc:
        exc.add_note(f"While calibrating {disp_file.name}")
        raise
    if result.unwrap_decisions is not None:
        decisions_csv = output_dir / "unwrap_decisions.csv"
        result.unwrap_decisions.to_csv(decisions_csv)
        logger.info(
            "Unwrap-error correction: %s (%s)",
            result.unwrap_decisions.summary(),
            decisions_csv,
        )
    _tropo_applied = result.components_applied["cal_tropo"]
    _set_applied = result.components_applied["cal_set"]
    cal_surface_m = result.calibration.astype(np.float32)
    nodata_pixel_count = int(np.isnan(cal_surface_m).sum())

    coords = {"time": time, "y": y, "x": x}
    calibration = xr.DataArray(
        cal_surface_m[np.newaxis, :, :],
        coords=coords,
        dims=["time", "y", "x"],
        attrs={
            "units": "meters",
            "long_name": "calibration_correction",
            "corrections_applied": (
                ", ".join(
                    name
                    for name, applied in [
                        ("troposphere", _tropo_applied),
                        ("solid_earth_tide", _set_applied),
                    ]
                    if applied
                )
                or "none"
            ),
        },
    )
    if result.sigma_cal is not None:
        # loclin surface: sigma_CAL (Venti R-E1) of the calibration itself
        std_values = np.asarray(result.sigma_cal, dtype=np.float32)
        s0_m = venti_options.uncertainty.sigma_nonsecular_meters
        uncertainty_source = (
            "sigma_CAL: sqrt((k * sigma_GNSS)^2 + sigma_fit^2"
            + (
                f" + sigma_nonsecular^2) with sigma_nonsecular = {s0_m * 1000:g} mm"
                if s0_m > 0
                else ")"
            )
            + f", k = {venti_options.uncertainty.k_grid} for the UNR"
            f" {cal.grid_type} grid and sigma_fit the kernel's standard error,"
            " grown inside interpolated defo/event areas"
            + (
                " and, where the fit had too little weight, extrapolated from the"
                " nearest fitted value and grown with the distance"
                if venti_options.uncertainty.extrapolate_sigma
                else ""
            )
            + ". It is the uncertainty"
            " of the calibration surface (the secular long-wavelength field)."
            " Add the DISP noise below the cutoff"
            + (
                "."
                if s0_m > 0
                else (
                    " and, for one pair, non-secular station motion (seasonal"
                    " loading; 2.6-4.9 mm on the validation frames)."
                )
            )
        )
    else:
        # windowed_plane (gamma): GNSS reference uncertainty, not that of the
        # fitted surface. Clipped at 0: Venti's RBF interpolation can overshoot.
        std_values = np.clip(gnss_los_std, 0, None)
        uncertainty_source = (
            "GNSS reference uncertainty: UNR"
            f" {cal.grid_type}-grid sigmas projected to LOS and interpolated"
            " to each pixel (rate sigma x interval for 'constant', combined"
            " reference/secondary position sigma for 'variable'). Not the"
            " uncertainty of the fitted calibration surface."
        )
    calibration_std = xr.DataArray(
        std_values[np.newaxis, :, :],
        coords=coords,
        dims=["time", "y", "x"],
        attrs={
            "units": "meters",
            "long_name": "calibration_uncertainty",
            "uncertainty_source": uncertainty_source,
        },
    )

    # /auxiliary: GNSS E/N/U displacement over the pair (zeros where the grid
    # type has no Venti sampling path, as before)
    coarse_y = y[::AUXILIARY_STEP]
    coarse_x = x[::AUXILIARY_STEP]
    coarse_coords = {"time": time, "y": coarse_y, "x": coarse_x}

    def _coarse_da(values: np.ndarray | None) -> xr.DataArray:
        if values is None:
            values = np.zeros((len(coarse_y), len(coarse_x)), dtype=np.float32)
        return xr.DataArray(
            values[np.newaxis, :, :], coords=coarse_coords, dims=["time", "y", "x"]
        )

    aux_comps = ("north_south", "east_west", "up_down")
    model_3d_da = {k: _coarse_da(model_3d[k] if model_3d else None) for k in aux_comps}
    model_3d_std_da = (
        {f"{k}_std": _coarse_da(model_3d_std[f"{k}_std"]) for k in aux_comps}
        if model_3d_std
        else None
    )
    # /plate_motion: rigid plate rotation over the pair, on the DISP grid.
    # Computed here, after the fit, from the LOS file in row windows, so it
    # adds nothing to the run's peak memory.
    plate = v_opts.frame.plate
    plate_motion = xr.DataArray(
        _plate_motion_los(
            plate,
            disp_product.get_epsg(),
            los_file,
            x,
            y,
            sec_decimal - ref_decimal,
        )[np.newaxis, :, :],
        coords=coords,
        dims=["time", "y", "x"],
    )
    logger.info("plate_motion: %s plate (ITRF2020-PMM)", plate)

    def _pkg_version(name: str) -> str:
        try:
            return importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            return "unknown"

    # The package's own version, not the installed distribution's: for an
    # editable install the dist-info is only refreshed on reinstall.
    cal_disp_version = __version__
    venti_version = _pkg_version("venti")

    _buf = StringIO()
    algorithm_parameters.to_yaml(_buf, with_comments=False)
    algorithm_parameters_yaml = _buf.getvalue()

    # Build CalProduct
    logger.info("Writing CalProduct to %s/...", output_dir)
    cal_product = CalProduct.create(
        calibration=calibration,
        disp_product=disp_product,
        output_dir=output_dir,
        calibration_std=calibration_std,
        spatial_ref=spatial_ref,
        global_metadata={
            # GNSS field is zero at the reference date
            "gnss_reference_epoch": f"{ref_decimal:.4f}",
            "auxiliary_model_3d_resolution": "5km",
            "calibration_resolution": f"{int(x_spacing)}m",
        },
        version=product_version,
        compression=compression,
    )

    cal_product.add_identification(
        calibration_reference_name=calibration_reference_name,
        calibration_reference_version=calibration_reference_version,
        calibration_reference_type=calibration_reference_type,
        calibration_reference_reference_frame=calibration_reference_reference_frame,
        source_data_file_list=source_data_file_list,
        source_calibration_file_list=source_calibration_file_list,
        source_data_access=source_data_access,
        source_data_dem_name=str(_meta("source_data_dem_name", "unknown")),
        source_data_satellite_names=source_data_satellite_names,
        source_data_imaging_geometry=str(
            _meta("source_data_imaging_geometry", "unknown")
        ),
        source_data_x_spacing=x_spacing,
        source_data_y_spacing=y_spacing,
        static_layers_data_access=static_layers_data_access,
        absolute_orbit_number=int(_meta("absolute_orbit_number", -1)),
        track_number=int(_meta("track_number", -1)),
        instrument_name=str(_meta("instrument_name", "unknown")),
        look_direction=str(_meta("look_direction", "unknown")).lower(),
        radar_band=str(_meta("radar_band", "unknown")),
        orbit_pass_direction=str(_meta("orbit_pass_direction", "unknown")).lower(),
        bounding_polygon=bounding_polygon,
        product_bounding_box=product_bounding_box,
        product_sample_spacing=product_sample_spacing,
        product_data_access=product_data_access,
        processing_facility=processing_facility,
        nodata_pixel_count=nodata_pixel_count,
        ceos_number_of_input_granules=len(source_data_file_list),
        processing_start_datetime=processing_start,
        ceos_analysis_ready_data_document_identifier=disp_meta.get(
            "ceos_analysis_ready_data_document_identifier"
        ),
        ceos_analysis_ready_data_product_type=disp_meta.get(
            "ceos_analysis_ready_data_product_type"
        ),
    )

    cal_product.add_metadata(
        algorithm_parameters_yaml=algorithm_parameters_yaml,
        platform_id=platform_id,
        source_data_software_disp_version=disp_product.version,
        cal_disp_software_version=cal_disp_version,
        venti_software_version=venti_version,
        product_pixel_coordinate_convention="center",
        ceos_atmospheric_phase_correction="tropospheric" if _tropo_applied else "none",
        ceos_gridding_convention="consistent",
        ceos_product_measurement_projection="line_of_sight",
        ceos_ionospheric_phase_correction="none",
        ceos_noise_removal="N",
        pge_runconfig=pge_runconfig,
        gnss_provenance_json=(
            json.dumps(gnss_provenance.as_dict())
            if gnss_provenance is not None
            else None
        ),
        # v0.5 route only, like gnss_provenance: the gamma product keeps its
        # 0.3 structure (T37.3)
        **(
            {
                "algorithm_schema_version": ALGORITHM_SCHEMA_VERSION,
                "frame_parameters_version": frame_parameters_version,
                "components_applied_json": json.dumps(result.components_applied),
            }
            if venti_options.surface.method == "loclin"
            else {}
        ),
    )

    span_years = sec_decimal - ref_decimal
    cal_product.add_auxiliary(
        model_3d=model_3d_da,
        model_3d_std=model_3d_std_da,
        spatial_ref=spatial_ref,
        pair_span_years=span_years,
        gnss_sampled=model_3d is not None,
    )
    cal_product.add_plate_motion(
        plate_motion,
        plate=plate,
        pair_span_years=span_years,
        spatial_ref=spatial_ref,
    )

    logger.info("Calibration complete: %s", cal_product.path)
    return cal_product.path
