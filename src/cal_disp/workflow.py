from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import xarray as xr

from cal_disp.product import CalProduct, DispProduct


def run_calibration(
    disp_file: Path,
    unr_grid_latlon_file: Path,  # noqa: ARG001 - TODO: add GNSS
    unr_timeseries_dir: Path,  # noqa: ARG001 - TODO: add GNSS
    output_dir: Path,
    dem_file: Path | None = None,  # noqa: ARG001 - TODO: add DEM
    los_file: Path | None = None,  # noqa: ARG001 - TODO: add LOS
    block_shape: tuple[int, int] = (
        512,
        512,
    ),  # noqa: ARG001 - TODO: implement blocking
    n_workers: int = 4,  # noqa: ARG001 - TODO: implement parallelization
    threads_per_worker: int = 1,  # noqa: ARG001 - TODO: implement parallelization
    work_directory: Path | None = None,  # noqa: ARG001 - TODO: implement temp files
    pge_runconfig: str | None = None,
    # Calibration reference parameters
    calibration_reference_name: str = "UNR gridded data",
    calibration_reference_version: str = "0.2",
    calibration_reference_type: str = "constant",
    calibration_reference_reference_frame: str = "IGS20",
    # Product metadata parameters
    platform_id: str = "S1A",
    absolute_orbit_number: int = 0,  # TODO: Extract from DISP product
    track_number: int = 0,  # TODO: Extract from DISP product
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
    # Software versions
    cal_disp_version: str = "0.1",
    venti_version: str = "0.1",
    disp_version: str = "1.0",
) -> Path:
    """Run displacement calibration workflow.

    Currently returns zeros - full implementation pending.

    Parameters
    ----------
    disp_file : Path
        Input displacement file.
    unr_grid_latlon_file : Path
        GNSS calibration latlon lookup for reference grid.
    unr_timeseries_dir : Path
        GNSS calibration reference grid dir with tenv8 files.
    output_dir : Path
        Output directory for calibrated displacement file.
    dem_file : Path or None, optional
        Digital elevation model file.
    los_file : Path or None, optional
        Line-of-sight geometry file.
    block_shape : tuple[int, int], optional
        Processing block size. Default is (512, 512).
    n_workers : int, optional
        Number of parallel workers. Default is 4.
    threads_per_worker : int, optional
        Threads per worker. Default is 1.
    work_directory : Path or None, optional
        Working directory for intermediate files.
    pge_runconfig : str or None, optional
        PGE run configuration string.
    calibration_reference_name : str, optional
        Calibration reference name. Default is "UNR gridded data".
    calibration_reference_version : str, optional
        Calibration reference version. Default is "0.2".
    calibration_reference_type : str, optional
        Calibration reference type. Default is "constant".
    calibration_reference_reference_frame : str, optional
        Reference frame. Default is "IGS20".
    platform_id : str, optional
        Platform ID. Default is "S1A".
    absolute_orbit_number : int, optional
        Absolute orbit number. Default is 0 (to be extracted).
    track_number : int, optional
        Track number. Default is 0 (to be extracted).
    instrument_name : str, optional
        Instrument name. Default is "C-SAR".
    look_direction : str, optional
        Look direction. Default is "right".
    radar_band : str, optional
        Radar band. Default is "C".
    orbit_pass_direction : str, optional
        Orbit direction. Default is "ascending".
    processing_facility : str, optional
        Processing facility. Default is "JPL".
    source_data_access : str, optional
        Source data access URL.
    source_data_dem_name : str, optional
        DEM name. Default is "Copernicus DEM GLO-30".
    source_data_imaging_geometry : str, optional
        Imaging geometry. Default is "right_looking".
    static_layers_data_access : str, optional
        Static layers URL.
    product_data_access : str, optional
        Product data access URL.
    cal_disp_version : str, optional
        cal_disp version. Default is "0.1".
    venti_version : str, optional
        venti version. Default is "0.1".
    disp_version : str, optional
        DISP processor version. Default is "1.0".

    Returns
    -------
    Path
        Path to output calibrated displacement file.

    Notes
    -----
    TODO: Implement actual calibration algorithm
    TODO: Load and use calibration_grid
    TODO: Apply DEM/LOS corrections
    TODO: Implement 3D velocity model estimation
    TODO: Extract orbit numbers and track from DISP metadata

    """
    # Load displacement product
    disp_product = DispProduct.from_path(disp_file)
    ds_disp = disp_product.open_dataset()

    # TBD
    block_shape = block_shape

    # Extract dimensions
    time = ds_disp.time.values
    y = ds_disp.y.values
    x = ds_disp.x.values
    shape = (len(time), len(y), len(x))

    # Extract spatial reference
    spatial_ref = ds_disp["spatial_ref"] if "spatial_ref" in ds_disp else None

    # Compute spatial metadata
    x_min, x_max = float(x.min()), float(x.max())
    y_min, y_max = float(y.min()), float(y.max())
    x_spacing = float(np.abs(x[1] - x[0]))
    y_spacing = float(np.abs(y[1] - y[0]))

    # Build bounding box and polygon
    product_bounding_box = f"({x_min}, {y_min}, {x_max}, {y_max})"
    bounding_polygon = (
        f"POLYGON(({x_min} {y_min}, {x_max} {y_min}, "
        f"{x_max} {y_max}, {x_min} {y_max}, {x_min} {y_min}))"
    )
    product_sample_spacing = f"{x_spacing}m"

    # Count nodata pixels
    nodata_pixel_count = int(np.isnan(ds_disp.displacement.values).sum())

    # Build file lists
    source_data_file_list = [str(disp_file.name)]

    # Build calibration file list from UNR inputs
    source_calibration_file_list = [str(unr_grid_latlon_file.name)]
    if unr_timeseries_dir.exists():
        tenv8_files = list(unr_timeseries_dir.glob("*.tenv8"))
        source_calibration_file_list.extend(
            [f.name for f in tenv8_files[:10]]
        )  # Limit to first 10

    # Determine satellite names from platform_id
    source_data_satellite_names = [
        f"Sentinel-{platform_id[-2:]}"
    ]  # e.g., "S1A" -> "Sentinel-1A"

    # Build coordinates for 3D data (time, y, x)
    coords = {"time": time, "y": y, "x": x}

    # TODO: Load GNSS calibration grid and compute actual calibration
    # For now, copy displacement values as placeholder
    calibration = xr.DataArray(
        ds_disp.displacement.values[np.newaxis, :, :],
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

    # TODO: Estimate 3D velocity model from GNSS
    # Coarse grid: every 167th point (~5km at 30m resolution)
    # Note: model_3d is 2D (y, x) - represents velocity rates, not time series
    coarse_y = y[::167]
    coarse_x = x[::167]
    coarse_shape = (len(time), len(coarse_y), len(coarse_x))

    model_3d = {
        "north_south": xr.DataArray(
            np.zeros(coarse_shape, dtype=np.float32),
            coords={"time": time, "y": coarse_y, "x": coarse_x},
            dims=["time", "y", "x"],
            attrs={"units": "meters/year", "long_name": "north_south_velocity"},
        ),
        "east_west": xr.DataArray(
            np.zeros(coarse_shape, dtype=np.float32),
            coords={"time": time, "y": coarse_y, "x": coarse_x},
            dims=["time", "y", "x"],
            attrs={"units": "meters/year", "long_name": "east_west_velocity"},
        ),
        "up_down": xr.DataArray(
            np.zeros(coarse_shape, dtype=np.float32),
            coords={"time": time, "y": coarse_y, "x": coarse_x},
            dims=["time", "y", "x"],
            attrs={"units": "meters/year", "long_name": "up_down_velocity"},
        ),
    }

    # Create calibrated product (main group only)
    cal = CalProduct.create(
        calibration=calibration,
        disp_product=disp_product,
        output_dir=output_dir,
        calibration_std=calibration_std,
        spatial_ref=spatial_ref,
        global_metadata={
            "gnss_reference_epoch": "2020-01-01T00:00:00Z",
            "auxiliary_model_3d_resolution": "5km",
            "calibration_resolution": "30m",
        },
        version="0.1",
    )

    # Add identification group
    cal.add_identification(
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
        processing_start_datetime=datetime.utcnow(),
    )

    # Add metadata group
    cal.add_metadata(
        algorithm_parameters_yaml="# TODO: Add actual algorithm parameters",
        platform_id=platform_id,
        source_data_software_disp_version=disp_version,
        cal_disp_software_version=cal_disp_version,
        venti_software_version=venti_version,
        product_pixel_coordinate_convention="center",
        ceos_atmospheric_phase_correction="none",
        ceos_gridding_convention="consistent",
        ceos_product_measurement_projection="line_of_sight",
        ceos_ionospheric_phase_correction="none",
        ceos_noise_removal="N",
        pge_runconfig=pge_runconfig,
    )

    # Add auxiliary group with 3D velocity model
    cal.add_auxiliary(
        model_3d=model_3d,
        spatial_ref=spatial_ref,
    )

    return cal.path
