"""Helpers shared by the DISP-CAL product writer tests."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

from cal_disp.product import CalProduct, DispProduct


def write_product(
    disp_file: Path, output_dir: Path, compression: bool = True, seed: int = 0
) -> tuple[CalProduct, np.ndarray]:
    """Write a complete product (all groups) from synthetic arrays.

    Returns the product and the ``(1, ny, nx)`` float32 calibration array
    that was written (with a NaN block in the upper-left corner).
    """
    disp = DispProduct.from_path(disp_file)
    with xr.open_dataset(disp_file, engine="h5netcdf") as ds:
        x, y, time = ds.x.values, ds.y.values, ds.time.values
        spatial_ref = ds["spatial_ref"].load()
    rng = np.random.default_rng(seed)
    cal = rng.normal(0, 0.05, (1, len(y), len(x))).astype(np.float32)
    cal[0, :10, :10] = np.nan
    coords = {"time": time, "y": y, "x": x}
    calibration = xr.DataArray(cal, coords=coords, dims=["time", "y", "x"])
    calibration_std = xr.DataArray(
        np.abs(cal) / 10, coords=coords, dims=["time", "y", "x"]
    )
    product = CalProduct.create(
        calibration=calibration,
        disp_product=disp,
        output_dir=output_dir,
        calibration_std=calibration_std,
        spatial_ref=spatial_ref,
        compression=compression,
    )
    product.add_identification(
        calibration_reference_name="UNR gridded data",
        calibration_reference_version="0.3",
        calibration_reference_type="constant",
        calibration_reference_reference_frame="IGS20",
        source_data_file_list=[disp_file.name],
        source_calibration_file_list=["grid_latlon_lookup_v0.3.txt"],
        source_data_access="https://example.org/disp",
        source_data_dem_name="Copernicus GLO-30",
        source_data_satellite_names=["S1A"],
        source_data_imaging_geometry="Geocoded",
        source_data_x_spacing=100.0,
        source_data_y_spacing=100.0,
        static_layers_data_access="https://example.org/static",
        absolute_orbit_number=1,
        track_number=2,
        instrument_name="C-SAR",
        look_direction="right",
        radar_band="C",
        orbit_pass_direction="ascending",
        bounding_polygon="POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))",
        product_bounding_box="(0, 0, 1, 1)",
        product_sample_spacing="100m",
        product_data_access="https://example.org/cal",
        processing_facility="test",
        nodata_pixel_count=100,
        ceos_number_of_input_granules=1,
        processing_start_datetime=datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc),
    )
    product.add_metadata(
        algorithm_parameters_yaml="a: 1\n",
        platform_id="S1A",
        source_data_software_disp_version="1.0",
        cal_disp_software_version="0.0",
        venti_software_version="0.0",
        pge_runconfig="cal_disp_workflow:\n  x: 1\n",
    )
    coarse = {"time": time, "y": y[::50], "x": x[::50]}
    zeros = np.zeros((1, len(coarse["y"]), len(coarse["x"])), dtype=np.float32)
    product.add_auxiliary(
        model_3d={
            comp: xr.DataArray(zeros, coords=coarse, dims=["time", "y", "x"])
            for comp in ("north_south", "east_west", "up_down")
        },
        spatial_ref=spatial_ref,
    )
    return product, cal
