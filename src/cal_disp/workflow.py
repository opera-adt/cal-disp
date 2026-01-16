from __future__ import annotations

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
    block_shape: tuple[int, int] = (  # noqa: ARG001 - TODO: implement blocking
        512,
        512,
    ),  # noqa: ARG001 - TODO: implement blocking
    n_workers: int = 4,  # noqa: ARG001 - TODO: implement parallelization
    threads_per_worker: int = 1,  # noqa: ARG001 - TODO: implement parallelization
    work_directory: Path | None = None,  # noqa: ARG001 - TODO: implement temp files
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

    """
    # Load displacement product
    disp_product = DispProduct.from_path(disp_file)
    ds_disp = disp_product.open_dataset()

    # Extract spatial dimensions
    y = ds_disp.y.values
    x = ds_disp.x.values
    shape = (len(y), len(x))

    # Extract spatial reference (it's a data variable, not a coordinate)
    spatial_ref = ds_disp["spatial_ref"] if "spatial_ref" in ds_disp else None

    # Build coordinates for 2D data (y, x only)
    # Note: time info preserved in metadata via CalProduct, not as coordinate
    coords = {"y": y, "x": x}

    # TODO: Load GNSS calibration grid and compute actual calibration
    # For now, create zeros (2D: y, x)
    calibration = xr.DataArray(
        ds_disp.displacement.values,
        coords=coords,
        dims=["y", "x"],
        attrs={"units": "meters", "long_name": "calibration_correction"},
    )

    calibration_std = xr.DataArray(
        np.zeros(shape, dtype=np.float32),
        coords=coords,
        dims=["y", "x"],
        attrs={"units": "meters", "long_name": "calibration_uncertainty"},
    )

    # TODO: Estimate 3D velocity model from GNSS
    # Coarse grid: every 167th point (~5km at 30m resolution)
    coarse_y = y[::167]
    coarse_x = x[::167]
    coarse_shape = (len(coarse_y), len(coarse_x))

    model_3d = {
        "north_south": xr.DataArray(
            np.zeros(coarse_shape, dtype=np.float32),
            coords={"y": coarse_y, "x": coarse_x},
            dims=["y", "x"],
            attrs={"units": "meters", "long_name": "north_south_displacement"},
        ),
        "east_west": xr.DataArray(
            np.zeros(coarse_shape, dtype=np.float32),
            coords={"y": coarse_y, "x": coarse_x},
            dims=["y", "x"],
            attrs={"units": "meters", "long_name": "east_west_displacement"},
        ),
        "up_down": xr.DataArray(
            np.zeros(coarse_shape, dtype=np.float32),
            coords={"y": coarse_y, "x": coarse_x},
            dims=["y", "x"],
            attrs={"units": "meters", "long_name": "up_down_displacement"},
        ),
    }

    # Create calibrated product with spatial reference
    cal = CalProduct.create(
        calibration=calibration,
        disp_product=disp_product,
        output_dir=output_dir,
        calibration_std=calibration_std,
        model_3d=model_3d,
        spatial_ref=spatial_ref,
        metadata={
            "gnss_reference_epoch": "2020-01-01T00:00:00Z",
            "auxiliary_model_3d_resolution": "5km",
            "calibration_resolution": "30m",
        },
    )

    return cal.path
