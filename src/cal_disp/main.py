from __future__ import annotations

from pathlib import Path

from cal_disp._log import get_max_memory_usage, log_runtime
from cal_disp._version import __version__
from cal_disp.browse_image import make_browse_image_from_nc
from cal_disp.config.workflow import CalibrationWorkflow
from cal_disp.workflow import run_calibration


@log_runtime
def run(runconfig: CalibrationWorkflow, debug: bool = False) -> Path:
    """Run the displacement calibration workflow.

    Parameters
    ----------
    runconfig : CalibrationWorkflow
        Workflow configuration for the calibration.
    debug : bool, optional
        Enable debug logging. Default is False.

    Returns
    -------
    Path
        Path to the output calibrated displacement file.

    Raises
    ------
    SystemExit
        If configuration validation fails or input files are missing.

    """
    # Setup logging
    logger = runconfig.setup_logging(level="DEBUG" if debug else "INFO")

    print(runconfig.summary())

    # Validate configuration
    status = runconfig.validate_ready_to_run()
    if not status["ready"]:
        logger.error("Configuration validation failed:")
        for err in status["errors"]:
            logger.error(f"  - {err}")
        raise SystemExit(1)

    runconfig.create_directories()
    logger.debug(f"Work directory: {runconfig.work_directory}")

    # Check input files exist
    missing = runconfig.get_missing_files()
    if missing:
        logger.error(f"Missing input files: {', '.join(missing)}")
        raise SystemExit(1)

    # Run calibration
    output_file = run_calibration(
        disp_file=Path(runconfig.input_options.disp_file),
        unr_grid_latlon_file=Path(runconfig.input_options.unr_grid_latlon_file),
        unr_timeseries_dir=Path(runconfig.input_options.unr_timeseries_dir),
        output_dir=runconfig.output_directory,
        dem_file=Path(runconfig.dynamic_ancillary_options.dem_file),
        los_file=Path(runconfig.dynamic_ancillary_options.los_file),
        block_shape=runconfig.worker_settings.block_shape,
        n_workers=runconfig.worker_settings.n_workers,
        threads_per_worker=runconfig.worker_settings.threads_per_worker,
        work_directory=runconfig.work_directory,
        pge_runconfig=str(runconfig._to_yaml_obj()),
    )

    # Generate browse image
    logger.info(f"Output product: {output_file}")
    output_png = output_file.with_suffix(".png")
    make_browse_image_from_nc(output_png, output_file)
    logger.info(f"Browse image: {output_png}")

    # Log metadata
    logger.info(f"cal_disp version: {__version__}")

    max_mem = get_max_memory_usage(units="GB")
    logger.info(f"Maximum memory usage: {max_mem:.2f} GB")

    return output_file
