"""Calibrate one OPERA DISP-S1 product against GNSS, using Venti.

Load the product, compute the GNSS LOS field from pre-staged UNR files,
estimate the calibration surface (``venti.estimate_calibration_surface``)
and write a ``CalProduct``.
"""

from __future__ import annotations

import functools
import importlib.metadata
import inspect
import logging
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path

import numpy as np
import rasterio
import xarray as xr

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
    displacement = ds_disp["displacement"]
    out = np.empty(displacement.shape, dtype=np.float32)
    for row in range(0, displacement.shape[0], block_rows):
        out[row : row + block_rows] = displacement[row : row + block_rows].values
    return out


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


def _find_reference_point(
    ds_disp: xr.Dataset, mask: np.ndarray, work_dir: Path
) -> tuple[int, int]:
    """Choose the reference pixel ``(row, col)`` for one product.

    Venti's rule (``opera_utils`` ``find_reference_point``) applied to the
    product's temporal coherence over valid pixels.
    """
    from opera_utils.disp import rebase_reference

    if "temporal_coherence" in ds_disp:
        quality = ds_disp.temporal_coherence.values.astype(np.float32)
    else:
        quality = np.ones(mask.shape, dtype=np.float32)
    quality = np.where(mask & np.isfinite(quality), quality, 0.0).astype(np.float32)

    x, y = ds_disp.x.values, ds_disp.y.values
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    crs_wkt = (
        ds_disp["spatial_ref"].attrs.get("crs_wkt")
        if "spatial_ref" in ds_disp
        else None
    )
    quality_file = work_dir / "reference_quality.tif"
    with rasterio.open(
        quality_file,
        "w",
        driver="GTiff",
        height=quality.shape[0],
        width=quality.shape[1],
        count=1,
        dtype="float32",
        crs=crs_wkt,
        transform=rasterio.transform.from_origin(
            float(x[0]) - dx / 2, float(y[0]) - dy / 2, dx, -dy
        ),
    ) as dst:
        dst.write(quality, 1)
    try:
        row, col = rebase_reference.find_reference_point(quality_file)
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
    block_shape: tuple[int, int] = (512, 512),  # noqa: ARG001 — TODO: Dask blocking
    n_workers: int = 4,  # noqa: ARG001 — TODO: Dask parallelisation
    threads_per_worker: int = 1,
    work_directory: Path | None = None,
    pge_runconfig: str | None = None,
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
    from venti.io import read_netcdf_correction
    from venti.surface import SENTINEL1_WAVELENGTH_M, estimate_calibration_surface

    ref_decimal = _date_to_decimal_year(disp_product.reference_date)
    gnss_los, gnss_los_std = _gnss_los_fields(
        gnss_ref,
        (los_east, los_north, los_up),
        disp_file,
        x,
        y,
        cal.downsample_factor,
        gnss_dir,
        ref_decimal,
        _date_to_decimal_year(disp_product.secondary_date),
    )

    # Displacement (m): float64 on disk, read in row blocks straight into
    # float32 so the float64 array is never held or cached in `ds_disp`
    disp_2d = _load_displacement(ds_disp)

    # Valid pixels (both masks: 1 = valid)
    mask = ~np.isnan(disp_2d)
    if "recommended_mask" in ds_disp:
        mask &= ds_disp.recommended_mask.values.astype(bool)
    if "water_mask" in ds_disp:
        mask &= ds_disp.water_mask.values.astype(bool)

    # Event areas (True = valid) are filled from neighbours before the fit
    event_mask = None
    if defo_area_db_json is not None or event_db_json is not None:
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

    # Signals GNSS does not contain: removed before the fit, added back to
    # the surface
    corrections: list[np.ndarray] = []
    _tropo_applied = False
    use_tropo = bool(reference_tropo_files and secondary_tropo_files)
    if use_tropo and not cal.apply_tropo_correction:
        logger.info("apply_tropo_correction=false; ignoring the tropo files")
        use_tropo = False
    if use_tropo:
        if dem_file is None:
            raise ValueError(
                "dem_file is required when tropospheric correction files are provided"
            )
        from cal_disp.prep.tropo import prepare_troposphere_correction

        logger.info("Preparing tropospheric correction (ref + sec)...")
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
        corrections.append(sec_tropo - ref_tropo)
        _tropo_applied = True

    _set_applied = False
    if cal.apply_solid_earth_tide_correction:
        set_corr = read_netcdf_correction(disp_file, "solid_earth_tide")
        if set_corr is None:
            logger.warning(
                "No /corrections/solid_earth_tide in %s; calibrating without SET"
                " correction",
                disp_file.name,
            )
        else:
            corrections.append(set_corr.astype(np.float32))
            _set_applied = True

    weights = None
    if cal.downsample_factor > 1 and cal.downsample_weighted:
        if "temporal_coherence" not in ds_disp:
            raise ValueError(
                "downsample_weighted=True but temporal_coherence is not in"
                f" {disp_file.name}"
            )
        weights = ds_disp.temporal_coherence.values.astype(np.float32)

    ref_point = _find_reference_point(ds_disp, mask, work_directory)
    # Cycle length of an unwrapping error: λ/2 (two-way path), derived from
    # the product's radar wavelength so other sensors (e.g. NISAR) are
    # handled too.  Venti calls this argument `wavelength_m` but uses it as
    # the value region offsets are rounded to.
    wavelength_m = (
        _read_wavelength_m(disp_file)
        if cal.unwrap_error_correction
        else SENTINEL1_WAVELENGTH_M
    )
    unwrap_cycle_m = _unwrap_cycle_length_m(wavelength_m)

    logger.info(
        "Fitting calibration surface (window=%d px, downsample=%d, corrections:"
        " tropo=%s, SET=%s)...",
        cal.window_size_pixels,
        cal.downsample_factor,
        _tropo_applied,
        _set_applied,
    )
    try:
        result = estimate_calibration_surface(
            disp_2d,
            gnss_los,
            mask,
            ref_point,
            cal.window_size_pixels,
            corrections=corrections,
            gnss_los_std=gnss_los_std if cal.weight_fit_by_gnss_uncertainty else None,
            event_mask=event_mask,
            options=cal.to_venti(),
            wavelength_m=unwrap_cycle_m,
            downsample_factor=cal.downsample_factor,
            downsample_method=cal.downsample_method,
            downsample_weights=weights,
            n_jobs=threads_per_worker,
        )
    except ValueError as exc:
        exc.add_note(f"While calibrating {disp_file.name}")
        raise
    cal_surface_m = result.surface.astype(np.float32)
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
    # TODO: improve. GNSS reference uncertainty, not that of the fitted
    # surface. Clipped at 0: Venti's RBF interpolation can overshoot below 0.
    calibration_std = xr.DataArray(
        np.clip(gnss_los_std, 0, None)[np.newaxis, :, :],
        coords=coords,
        dims=["time", "y", "x"],
        attrs={
            "units": "meters",
            "long_name": "calibration_uncertainty",
            "uncertainty_source": (
                "GNSS reference uncertainty: UNR"
                f" {cal.grid_type}-grid sigmas projected to LOS and interpolated"
                " to each pixel (rate sigma x interval for 'constant', combined"
                " reference/secondary position sigma for 'variable'). Not the"
                " uncertainty of the fitted calibration surface."
            ),
        },
    )

    # 3-D velocity model: zero placeholder until decomposition exists
    coarse_y = y[::167]
    coarse_x = x[::167]
    coarse_shape = (len(time), len(coarse_y), len(coarse_x))
    coarse_coords = {"time": time, "y": coarse_y, "x": coarse_x}

    def _zero_da(long_name: str, units: str) -> xr.DataArray:
        return xr.DataArray(
            np.zeros(coarse_shape, dtype=np.float32),
            coords=coarse_coords,
            dims=["time", "y", "x"],
            attrs={"units": units, "long_name": long_name},
        )

    model_3d = {
        "north_south": _zero_da("north_south_velocity", "meters/year"),
        "east_west": _zero_da("east_west_velocity", "meters/year"),
        "up_down": _zero_da("up_down_velocity", "meters/year"),
    }

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
    )

    cal_product.add_auxiliary(
        model_3d=model_3d,
        spatial_ref=spatial_ref,
    )

    logger.info("Calibration complete: %s", cal_product.path)
    return cal_product.path
