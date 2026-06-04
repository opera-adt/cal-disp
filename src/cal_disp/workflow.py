"""Calibration workflow for OPERA DISP-S1 products.

Orchestrates GNSS-based calibration using Venti as the processing back-end:

1. Load the DISP product and extract spatial metadata.
2. Set up a ``GNSSReference`` pointing at pre-staged UNR tenv8 files.
3. Load the 3-band LOS (east / north / up) GeoTIFF.
4. Compute the GNSS LOS reference (velocity or epoch-displacement).
5. Fit a calibration surface with ``SpatialProcessor.fit_windowed_surface``.
6. Package the result into a ``CalProduct`` NetCDF.
"""

from __future__ import annotations

import importlib.metadata
import logging
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path

import numpy as np
import rasterio
import xarray as xr

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
    """Generate and combine event/deformation masks from GeoJSON database.

    For each available GeoJSON, selects only features whose ``frame_id``
    and ``event_date`` overlap the DISP product epoch, and ANDs the results
    into a single combined mask GeoTIFF (1 = valid, 0 = masked).
    Returns ``None`` when no database is provided.
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
    """Read the radar wavelength in meters from a DISP product file.

    Returns ~0.05546 for Sentinel-1 C-band.
    """
    from netCDF4 import Dataset  # type: ignore[import-untyped]

    try:
        with Dataset(disp_file, "r") as nc:
            ident_group = nc.groups["identification"]
            wl_m = float(ident_group.variables["radar_wavelength"][:])
            wl_units = getattr(ident_group.variables["radar_wavelength"], "units", "m")
        logger.debug(
            "Read radar_wavelength %.6f %s from %s", wl_m, wl_units, disp_file.name
        )
        return wl_m
    except Exception as e:
        msg = (
            f"Could not read /identification/radar_wavelength from {disp_file.name}. "
            "Ensure the file is a valid DISP product with an /identification group."
        )
        raise RuntimeError(msg) from e


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
    """Load a 3-band LOS GeoTIFF and return (east, north, up) arrays."""
    with rasterio.open(los_file) as src:
        if src.count < 3:
            raise ValueError(
                "LOS file must have ≥3 bands (east, north, up),"
                f" got {src.count}: {los_file}"
            )
        los_east = src.read(1).astype(np.float32)
        los_north = src.read(2).astype(np.float32)
        los_up = src.read(3).astype(np.float32)
    return los_east, los_north, los_up


def _setup_gnss_reference(
    disp_product: DispProduct,
    unr_grid_latlon_file: Path,
    unr_timeseries_dir: Path,
    gnss_dir: Path,
    reference_frame: str,
):
    """Initialise a ``GNSSReference`` reusing pre-staged UNR files.

    Cal-disp pre-stages UNR tenv8 files via ``cal-disp download unr``.
    Venti's ``GNSSReference.download_stations()`` will skip any file already
    present in ``output_dir``, so we symlink the staged files into a dedicated
    gnss working directory before calling it.
    """
    from venti.gnss import GNSSReference

    gnss_dir.mkdir(parents=True, exist_ok=True)

    # Venti expects the lookup at output_dir / "grid_latlon_lookup.txt"
    lookup_link = gnss_dir / "grid_latlon_lookup.txt"
    if not lookup_link.exists() and unr_grid_latlon_file.exists():
        lookup_link.symlink_to(unr_grid_latlon_file.resolve())

    # Symlink pre-staged tenv8 files so download_stations() uses them directly
    for tenv8 in unr_timeseries_dir.glob("*.tenv8"):
        link = gnss_dir / tenv8.name
        if not link.exists():
            link.symlink_to(tenv8.resolve())

    # Extract UTM bounds (south, north, west, east) from the DISP product
    b = disp_product.get_bounds()
    bounds = (b["bottom"], b["top"], b["left"], b["right"])
    utm_epsg = disp_product.get_epsg()

    gnss_ref = GNSSReference(
        bounds=bounds,
        output_dir=gnss_dir,
        reference_frame=reference_frame,
        utm_epsg=utm_epsg,
    )
    n_stations = gnss_ref.download_stations()
    logger.info("GNSS setup complete: %d stations available", n_stations)
    return gnss_ref


# Public entry point
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
    calibration_reference_type: str = "constant",
    calibration_reference_reference_frame: str = "IGS20",
    # Product metadata
    platform_id: str = "S1A",
    absolute_orbit_number: int = 0,  # TODO: extract from DISP identification group
    track_number: int = 0,  # TODO: extract from DISP identification group
    instrument_name: str = "C-SAR",
    look_direction: str = "right",
    radar_band: str = "C",
    orbit_pass_direction: str = "ascending",
    processing_facility: str = "JPL",
    source_data_access: str = "https://datapool.asf.alaska.edu/",
    source_data_dem_name: str = "Copernicus DEM GLO-30",
    source_data_imaging_geometry: str = "right_looking",
    static_layers_data_access: str = "https://example.com/static_layers",
    product_data_access: str = "https://example.com/products",
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
        DEM GeoTIFF.  Required when tropospheric correction files are provided.
    reference_tropo_files : list[Path], optional
        OPERA TROPO-ZENITH NetCDF files for the reference date (1–2 files for
        temporal interpolation).  When provided together with
        ``secondary_tropo_files``, a differential LOS correction
        (secondary − reference) is applied to the displacement before surface
        fitting.
    secondary_tropo_files : list[Path], optional
        OPERA TROPO-ZENITH NetCDF files for the secondary date (1–2 files for
        temporal interpolation).
    los_file : Path
        3-band LOS GeoTIFF (band 1 = east, 2 = north, 3 = up unit vectors).
    defo_area_db_json : Path, optional
        GeoJSON of continuous deformation areas to exclude from calibration.
        Features are filtered by ``frame_id`` and ``event_date`` matching the
        DISP product epoch.
    event_db_json : Path, optional
        GeoJSON of earthquake/volcanic events to exclude from calibration.
        Same filtering logic as ``defo_area_db_json``.
    block_shape : tuple[int, int]
        Dask block shape — reserved for future parallelisation.
    n_workers : int
        Dask workers — reserved for future parallelisation.
    threads_per_worker : int
        Threads per Dask worker, forwarded to ``SpatialProcessor`` ``n_jobs``.
    work_directory : Path, optional
        Scratch directory for intermediate files and GNSS cache.
    pge_runconfig : str, optional
        Serialised PGE run-config YAML stored in product metadata.
    calibration_reference_name : str
        Human-readable name of the calibration reference dataset.
    calibration_reference_version : str
        Version string of the calibration reference dataset.
    calibration_reference_type : str
        ``'constant'`` or ``'variable'``.
    calibration_reference_reference_frame : str
        GNSS reference frame (e.g. ``'IGS20'``).
    platform_id : str
        Satellite platform identifier (e.g. ``'S1A'``).
    absolute_orbit_number : int
        Absolute orbit number (extracted from DISP product in future).
    track_number : int
        Track number (extracted from DISP product in future).
    instrument_name : str
        SAR instrument name.
    look_direction : str
        Radar look direction.
    radar_band : str
        Radar frequency band.
    orbit_pass_direction : str
        Orbit pass direction (``'ascending'`` or ``'descending'``).
    processing_facility : str
        Processing facility name.
    source_data_access : str
        URL for source data access.
    source_data_dem_name : str
        Name of the DEM used.
    source_data_imaging_geometry : str
        Imaging geometry description.
    static_layers_data_access : str
        URL for static layer access.
    product_data_access : str
        URL for product data access.

    Returns
    -------
    Path
        Path to the output ``CalProduct`` NetCDF file.

    """
    from venti.spatial import SpatialProcessor, downsample_array, upsample_array

    if los_file is None:
        raise ValueError("los_file is required for calibration")

    if algorithm_parameters is None:
        algorithm_parameters = AlgorithmParameters()

    if work_directory is None:
        work_directory = output_dir / "scratch"
    work_directory.mkdir(parents=True, exist_ok=True)

    cal = algorithm_parameters.calibration_options

    # Load DISP product
    logger.info("Loading DISP product: %s", disp_file.name)
    disp_product = DispProduct.from_path(disp_file)
    ds_disp = disp_product.open_dataset()

    time = ds_disp.time.values
    y = ds_disp.y.values
    x = ds_disp.x.values
    shape = (len(time), len(y), len(x))

    spatial_ref = ds_disp.get("spatial_ref")

    x_min, x_max = float(x.min()), float(x.max())
    y_min, y_max = float(y.min()), float(y.max())
    x_spacing = float(np.abs(x[1] - x[0]))
    y_spacing = float(np.abs(y[1] - y[0]))

    product_bounding_box = f"({x_min}, {y_min}, {x_max}, {y_max})"
    bounding_polygon = (
        f"POLYGON(({x_min} {y_min}, {x_max} {y_min}, "
        f"{x_max} {y_max}, {x_min} {y_max}, {x_min} {y_min}))"
    )
    product_sample_spacing = f"{x_spacing}m"
    nodata_pixel_count = int(np.isnan(ds_disp.displacement.values).sum())

    source_data_file_list = [disp_file.name]
    source_calibration_file_list = [unr_grid_latlon_file.name]
    tenv8_files = sorted(unr_timeseries_dir.glob("*.tenv8"))
    source_calibration_file_list.extend(f.name for f in tenv8_files[:10])

    source_data_satellite_names = [f"Sentinel-{platform_id[-2:]}"]

    # GNSS reference setup
    logger.info("Setting up GNSS reference...")
    gnss_dir = work_directory / "gnss"
    gnss_ref = _setup_gnss_reference(
        disp_product=disp_product,
        unr_grid_latlon_file=unr_grid_latlon_file,
        unr_timeseries_dir=unr_timeseries_dir,
        gnss_dir=gnss_dir,
        reference_frame=cal.reference_frame,
    )

    # Load LOS unit vectors
    logger.info("Loading LOS unit vectors: %s", los_file.name)
    los_east, los_north, los_up = _load_los_bands(los_file)

    # Compute GNSS LOS reference
    from venti.gnss import compute_gnss_los

    ref_decimal = _date_to_decimal_year(disp_product.reference_date)
    sec_decimal = _date_to_decimal_year(disp_product.secondary_date)

    gnss_los = (
        compute_gnss_los(
            gnss_ref=gnss_ref,
            los_east=los_east,
            los_north=los_north,
            los_up=los_up,
            netcdf_file=disp_file,
            grid_type=cal.grid_type,
            cache_dir=gnss_dir,
            ref_date=ref_decimal,
            sec_date=sec_decimal,
            starting_year=cal.starting_year,
        )
        / 1000.0
    )

    # Prepare displacement array.
    # A single DISP file encodes one (ref_date, sec_date) pair.
    # displacement is 2-D (y, x); time is a length-1 coordinate.
    disp_2d = ds_disp.displacement.values.astype(np.float32)

    # Build valid-pixel mask
    mask = ~np.isnan(disp_2d)
    if "recommended_mask" in ds_disp:
        # recommended_mask convention: 1 = valid/recommended, 0 = invalid
        mask &= ds_disp.recommended_mask.values.astype(bool)
    if "water_mask" in ds_disp:
        # water_mask convention: 1 = land/valid, 0 = water/invalid
        mask &= ds_disp.water_mask.values.astype(bool)
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
            mask &= event_mask
            logger.info(
                "Applied event mask: %d pixels excluded", int((~event_mask).sum())
            )

    disp_masked = np.where(mask, disp_2d, np.nan)

    # Optional unwrap-error correction
    if cal.unwrap_error_correction:
        from venti.unwrap import correct_region_offset

        wavelength_m = _read_wavelength_m(disp_file)
        logger.info(
            "Applying unwrap-error correction (wavelength=%.5f m)...",
            wavelength_m,
        )
        disp_masked = correct_region_offset(
            input_disp=disp_masked,
            mask=mask,
            wavelength=wavelength_m,
        )
        if isinstance(disp_masked, np.ma.MaskedArray):
            disp_masked = disp_masked.filled(np.nan)

    # Fit calibration surface
    spatial_processor = SpatialProcessor()
    win_px = cal.window_size_pixels
    ds_factor = cal.downsample_factor
    original_shape = disp_masked.shape

    if ds_factor > 1:
        weights = None
        if cal.downsample_weighted:
            if "temporal_coherence" in ds_disp:
                weights = ds_disp.temporal_coherence.values.astype(np.float32)
                logger.info(
                    "Downsampling weighted by temporal coherence (factor=%d)", ds_factor
                )
            else:
                logger.warning(
                    "downsample_weighted=True but temporal_coherence not found in DISP "
                    "product; falling back to unweighted downsampling"
                )
        else:
            logger.info(
                "Downsampling by factor %d (%s)...", ds_factor, cal.downsample_method
            )

        disp_ds = downsample_array(
            disp_masked, ds_factor, method=cal.downsample_method, weights=weights
        )
        gnss_ds = downsample_array(gnss_los, ds_factor, method=cal.downsample_method)
        win_x = max(1, win_px // ds_factor)
        win_y = max(1, win_px // ds_factor)
    else:
        disp_ds = disp_masked
        gnss_ds = gnss_los
        win_x = win_px
        win_y = win_px

    logger.info(
        "Fitting calibration surface (window=%d×%d px, overlap=50%%, smoothing=%s)...",
        win_x,
        win_y,
        cal.calibration_surface_smoothing_method,
    )
    cal_surface = spatial_processor.fit_windowed_surface(
        insar_data=disp_ds,
        gnss_los=gnss_ds,
        window_size_x=win_x,
        window_size_y=win_y,
        window_overlap_x=win_x // 2,
        window_overlap_y=win_y // 2,
        n_jobs=threads_per_worker,
        smoothing_sigma=cal.calibration_surface_smoothing_sigma,
        smoothing_method=cal.calibration_surface_smoothing_method,
        sg_window_length=cal.savitzky_golay.window_length,
        sg_polyorder=cal.savitzky_golay.polyorder,
    )

    if ds_factor > 1:
        cal_surface = upsample_array(cal_surface, original_shape)

    cal_surface_m = cal_surface.astype(np.float32)

    # Optional tropospheric correction: prepare and apply to calibration surface
    _tropo_applied = False
    tropo_corr: np.ndarray | None = None
    if reference_tropo_files and secondary_tropo_files:
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
        tropo_corr = sec_tropo - ref_tropo  # differential LOS delay (meters)
        cal_surface_m = cal_surface_m + tropo_corr
        _tropo_applied = True
        logger.info("Tropospheric correction applied to calibration surface.")

    # Insert time dimension

    coords = {"time": time, "y": y, "x": x}
    calibration = xr.DataArray(
        cal_surface_m[np.newaxis, :, :],
        coords=coords,
        dims=["time", "y", "x"],
        attrs={"units": "meters", "long_name": "calibration_correction"},
    )
    calibration_std = xr.DataArray(
        np.zeros(shape, dtype=np.float32),
        coords=coords,
        dims=["time", "y", "x"],
        attrs={"units": "meters", "long_name": "calibration_uncertainty"},
    )

    # Coarse 3-D velocity model (placeholder — DecompositionWorkflow TBD)
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

    # Software version strings
    def _pkg_version(name: str) -> str:
        try:
            return importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            return "unknown"

    cal_disp_version = _pkg_version("cal-disp")
    venti_version = _pkg_version("venti")

    # Serialise algorithm parameters to YAML string for metadata
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
            "gnss_reference_epoch": f"{cal.starting_year:.1f}",
            "auxiliary_model_3d_resolution": "5km",
            "calibration_resolution": f"{int(x_spacing)}m",
        },
        version="0.1",
    )

    cal_product.add_identification(
        calibration_reference_name=calibration_reference_name,
        calibration_reference_version=calibration_reference_version,
        calibration_reference_type=calibration_reference_type,
        calibration_reference_reference_frame=calibration_reference_reference_frame,
        source_data_file_list=source_data_file_list,
        source_calibration_file_list=source_calibration_file_list,
        source_data_access=source_data_access,
        source_data_dem_name=source_data_dem_name,
        source_data_satellite_names=source_data_satellite_names,
        source_data_imaging_geometry=source_data_imaging_geometry,
        source_data_x_spacing=x_spacing,
        source_data_y_spacing=y_spacing,
        static_layers_data_access=static_layers_data_access,
        absolute_orbit_number=absolute_orbit_number,
        track_number=track_number,
        instrument_name=instrument_name,
        look_direction=look_direction,
        radar_band=radar_band,
        orbit_pass_direction=orbit_pass_direction,
        bounding_polygon=bounding_polygon,
        product_bounding_box=product_bounding_box,
        product_sample_spacing=product_sample_spacing,
        product_data_access=product_data_access,
        processing_facility=processing_facility,
        nodata_pixel_count=nodata_pixel_count,
        ceos_number_of_input_granules=len(source_data_file_list),
        processing_start_datetime=datetime.now(tz=timezone.utc),
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
