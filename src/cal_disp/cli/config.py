from __future__ import annotations

import functools
from pathlib import Path
from typing import Final

import click

# Show defaults for all options
click.option = functools.partial(click.option, show_default=True)

# Defaults
DEFAULT_CONFIG_NAME: Final[str] = "runconfig.yaml"
DEFAULT_N_WORKERS: Final[int] = 4
DEFAULT_THREADS_PER_WORKER: Final[int] = 1
DEFAULT_BLOCK_SHAPE: Final[tuple[int, int]] = (512, 512)
DEFAULT_PRODUCT_VERSION: Final[str] = "1.0"


def create_config(
    disp_file: Path,
    calibration_grid_latlon_file: Path,
    calibration_grid_ts_dir: Path,
    frame_id: int,
    algorithm_params_file: Path,
    los_file: Path,
    dem_file: Path,
    output_dir: Path,
    work_dir: Path = Path.cwd(),
    config_name: str = DEFAULT_CONFIG_NAME,
    mask_file: Path | None = None,
    ref_tropo_files: list[Path] | None = None,
    sec_tropo_files: list[Path] | None = None,
    iono_files: list[Path] | None = None,
    tiles_files: list[Path] | None = None,
    algorithm_overrides_json: Path | None = None,
    defo_area_db_json: Path | None = None,
    event_db_json: Path | None = None,
    n_workers: int = DEFAULT_N_WORKERS,
    threads_per_worker: int = DEFAULT_THREADS_PER_WORKER,
    block_shape: tuple[int, int] = DEFAULT_BLOCK_SHAPE,
    product_version: str = DEFAULT_PRODUCT_VERSION,
    output_format: str = "netcdf",
    compression: bool = True,
    keep_relative: bool = False,
) -> Path:
    """Generate and save a run configuration file for displacement calibration.

    Parameters
    ----------
    disp_file : Path
        Input displacement file to calibrate.
    calibration_grid_latlon_file : Path
        UNR grid lookup table (grid_latlon_lookup_v0.2.txt).
    calibration_grid_ts_dir : Path
        Directory containing UNR .tenv8 timeseries files.
    frame_id : int
        Frame ID of the DISP frame.
    algorithm_params_file : Path
        Algorithm parameters configuration file.
    los_file : Path
        DISP static LOS layer file (line-of-sight unit vectors).
    dem_file : Path
        DISP static DEM layer file (digital elevation model).
    output_dir : Path
        Directory for output files.
    work_dir : Path, optional
        Working directory for temporary files. Default is current directory.
    config_name : str, optional
        Name of configuration file to create. Default is "runconfig.yaml".
    mask_file : Path or None, optional
        Byte mask file to ignore low correlation/bad data (0=invalid, 1=good).
    ref_tropo_files : list[Path] or None, optional
        TROPO files for reference (primary) date.
    sec_tropo_files : list[Path] or None, optional
        TROPO files for secondary date.
    iono_files : list[Path] or None, optional
        Ionospheric correction files.
    tiles_files : list[Path] or None, optional
        Calibration tile bounds files (e.g., S1 burst bounds).
    algorithm_overrides_json : Path or None, optional
        Frame-specific algorithm parameter overrides.
    defo_area_db_json : Path or None, optional
        GeoJSON with deforming areas to exclude from calibration.
    event_db_json : Path or None, optional
        GeoJSON with earthquake/volcanic activity events.
    n_workers : int, optional
        Number of parallel workers. Default is 4.
    threads_per_worker : int, optional
        Threads per worker. Default is 1.
    block_shape : tuple[int, int], optional
        Processing block size (height, width). Default is (512, 512).
    product_version : str, optional
        Output product version. Default is "1.0".
    output_format : str, optional
        Output file format. Default is "netcdf".
    compression : bool, optional
        Whether to compress output. Default is True.
    keep_relative : bool, optional
        Keep paths relative instead of absolute. Default is False.

    Returns
    -------
    Path
        Path to created configuration file.

    Raises
    ------
    FileNotFoundError
        If required input files don't exist.
    ValueError
        If configuration parameters are invalid.

    """
    from cal_disp.config import (
        DynamicAncillaryFileGroup,
        InputFileGroup,
        StaticAncillaryFileGroup,
        WorkerSettings,
    )
    from cal_disp.config.pge_runconfig import (
        OutputOptions,
        PrimaryExecutable,
        ProductPathGroup,
        RunConfig,
    )

    # Validate required inputs exist
    required_files = {
        "disp_file": disp_file,
        "calibration_grid_latlon_file": calibration_grid_latlon_file,
        "algorithm_params_file": algorithm_params_file,
        "los_file": los_file,
        "dem_file": dem_file,
    }
    for name, path in required_files.items():
        if not path.exists():
            raise FileNotFoundError(f"{name} not found: {path}")

    # Validate calibration grid directory
    if not calibration_grid_ts_dir.exists():
        raise FileNotFoundError(
            f"calibration_grid_ts_dir not found: {calibration_grid_ts_dir}"
        )
    if not calibration_grid_ts_dir.is_dir():
        raise ValueError(
            f"calibration_grid_ts_dir must be a directory: {calibration_grid_ts_dir}"
        )

    # Validate optional single files if provided
    optional_files = {
        "mask_file": mask_file,
        "algorithm_overrides_json": algorithm_overrides_json,
        "defo_area_db_json": defo_area_db_json,
        "event_db_json": event_db_json,
    }
    for name, path in optional_files.items():  # type: ignore[assignment]
        if path is None:
            continue
        if not path.exists():
            raise FileNotFoundError(f"{name} not found: {path}")

    # Validate optional file lists if provided
    for name, file_list in [
        ("ref_tropo_files", ref_tropo_files),
        ("sec_tropo_files", sec_tropo_files),
        ("iono_files", iono_files),
        ("tiles_files", tiles_files),
    ]:
        if file_list:
            for f in file_list:
                if not f.exists():
                    raise FileNotFoundError(f"{name} file not found: {f}")

    # Validate parameters
    if frame_id < 0:
        raise ValueError(f"frame_id must be >= 0, got {frame_id}")
    if n_workers < 1:
        raise ValueError(f"n_workers must be >= 1, got {n_workers}")
    if threads_per_worker < 1:
        raise ValueError(f"threads_per_worker must be >= 1, got {threads_per_worker}")
    if len(block_shape) != 2 or any(dim < 1 for dim in block_shape):
        raise ValueError(
            f"block_shape must be (height, width) with positive ints, got {block_shape}"
        )

    # Create directories
    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    config_path = work_dir / config_name

    # Resolve to absolute paths unless keeping relative
    if not keep_relative:
        disp_file = disp_file.resolve()
        calibration_grid_latlon_file = calibration_grid_latlon_file.resolve()
        calibration_grid_ts_dir = calibration_grid_ts_dir.resolve()
        algorithm_params_file = algorithm_params_file.resolve()
        los_file = los_file.resolve()
        dem_file = dem_file.resolve()
        output_dir = output_dir.resolve()
        work_dir = work_dir.resolve()

        if mask_file is not None:
            mask_file = mask_file.resolve()
        if algorithm_overrides_json is not None:
            algorithm_overrides_json = algorithm_overrides_json.resolve()
        if defo_area_db_json is not None:
            defo_area_db_json = defo_area_db_json.resolve()
        if event_db_json is not None:
            event_db_json = event_db_json.resolve()

        # Resolve file lists
        if ref_tropo_files is not None:
            ref_tropo_files = [f.resolve() for f in ref_tropo_files]
        if sec_tropo_files is not None:
            sec_tropo_files = [f.resolve() for f in sec_tropo_files]
        if iono_files is not None:
            iono_files = [f.resolve() for f in iono_files]
        if tiles_files is not None:
            tiles_files = [f.resolve() for f in tiles_files]

    # Build configuration groups
    input_file_group = InputFileGroup(
        disp_file=disp_file,
        calibration_reference_latlon_file=calibration_grid_latlon_file,
        calibration_reference_grid_dir=calibration_grid_ts_dir,
        frame_id=frame_id,
    )

    dynamic_ancillary_group = DynamicAncillaryFileGroup(
        algorithm_parameters_file=algorithm_params_file,
        los_file=los_file,
        dem_file=dem_file,
        mask_file=mask_file,
        reference_tropo_files=ref_tropo_files,
        secondary_tropo_files=sec_tropo_files,
        iono_files=iono_files,
        tiles_files=tiles_files,
    )

    # Only create static ancillary group if any files are provided
    static_ancillary_group = None
    if any([algorithm_overrides_json, defo_area_db_json, event_db_json]):
        static_ancillary_group = StaticAncillaryFileGroup(
            algorithm_parameters_overrides_json=algorithm_overrides_json,
            deformation_area_database_json=defo_area_db_json,
            event_database_json=event_db_json,
        )

    product_path_group = ProductPathGroup(
        product_path=output_dir,
        scratch_path=work_dir,
        output_path=output_dir,
    )

    worker_settings = WorkerSettings(
        n_workers=n_workers,
        threads_per_worker=threads_per_worker,
        block_shape=block_shape,
    )

    output_options = OutputOptions(
        product_version=product_version,
        output_format=output_format,
        compression=compression,
    )

    # Create runconfig
    runconfig = RunConfig(
        input_file_group=input_file_group,
        dynamic_ancillary_group=dynamic_ancillary_group,
        static_ancillary_group=static_ancillary_group,
        output_options=output_options,
        primary_executable=PrimaryExecutable(),
        product_path_group=product_path_group,
        worker_settings=worker_settings,
    )

    runconfig.to_yaml(config_path)
    return config_path


@click.command("config")
@click.option(
    "--config-file",
    "-c",
    type=click.Path(path_type=Path),
    default=Path.cwd() / DEFAULT_CONFIG_NAME,
    help="Output configuration file path.",
)
@click.option(
    "--disp-file",
    "-d",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Input displacement file.",
)
@click.option(
    "--calibration-grid-latlon",
    "-cl",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="UNR grid lookup table (grid_latlon_lookup_v0.2.txt).",
)
@click.option(
    "--calibration-grid-dir",
    "-cd",
    type=click.Path(exists=True, path_type=Path, file_okay=False, dir_okay=True),
    required=True,
    help="Directory containing UNR .tenv8 timeseries files.",
)
@click.option(
    "--frame-id",
    "-f",
    type=int,
    required=True,
    help="Frame ID of the DISP frame.",
)
@click.option(
    "--algorithm-params",
    "-a",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="Algorithm parameters configuration file.",
)
@click.option(
    "--los-file",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="DISP static LOS layer file.",
)
@click.option(
    "--dem-file",
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help="DISP static DEM layer file.",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(path_type=Path),
    required=True,
    help="Output directory for calibrated products.",
)
@click.option(
    "--work-dir",
    type=click.Path(path_type=Path),
    default=Path.cwd(),
    help="Working directory for temporary files and scratch data.",
)
@click.option(
    "--mask-file",
    type=click.Path(exists=True, path_type=Path),
    help="Byte mask file (0=invalid, 1=good).",
)
@click.option(
    "--ref-tropo-files",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="TROPO files for reference date (can specify multiple times).",
)
@click.option(
    "--sec-tropo-files",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="TROPO files for secondary date (can specify multiple times).",
)
@click.option(
    "--iono-files",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="Ionospheric correction files (can specify multiple times).",
)
@click.option(
    "--tiles-files",
    type=click.Path(exists=True, path_type=Path),
    multiple=True,
    help="Calibration tile bounds files (can specify multiple times).",
)
@click.option(
    "--algorithm-overrides",
    type=click.Path(exists=True, path_type=Path),
    help="Frame-specific algorithm parameter overrides (JSON).",
)
@click.option(
    "--defo-area-db",
    type=click.Path(exists=True, path_type=Path),
    help="Deforming areas database (GeoJSON).",
)
@click.option(
    "--event-db",
    type=click.Path(exists=True, path_type=Path),
    help="Events database (GeoJSON).",
)
@click.option(
    "--n-workers",
    "-w",
    type=int,
    default=DEFAULT_N_WORKERS,
    help="Number of parallel workers.",
)
@click.option(
    "--threads-per-worker",
    "-t",
    type=int,
    default=DEFAULT_THREADS_PER_WORKER,
    help="Threads per worker.",
)
@click.option(
    "--block-shape",
    type=(int, int),
    default=DEFAULT_BLOCK_SHAPE,
    help="Processing block shape (height, width).",
)
@click.option(
    "--product-version",
    type=str,
    default=DEFAULT_PRODUCT_VERSION,
    help="Output product version.",
)
@click.option(
    "--output-format",
    type=click.Choice(["netcdf", "hdf5"]),
    default="netcdf",
    help="Output file format.",
)
@click.option(
    "--compression/--no-compression",
    default=True,
    help="Enable/disable output compression.",
)
@click.option(
    "--keep-relative/--absolute-paths",
    default=False,
    help="Keep paths relative instead of absolute.",
)
def config_cli(
    config_file: Path,
    disp_file: Path,
    calibration_grid_latlon: Path,
    calibration_grid_dir: Path,
    frame_id: int,
    algorithm_params: Path,
    los_file: Path,
    dem_file: Path,
    output_dir: Path,
    work_dir: Path,
    mask_file: Path | None,
    ref_tropo_files: tuple[Path, ...],
    sec_tropo_files: tuple[Path, ...],
    iono_files: tuple[Path, ...],
    tiles_files: tuple[Path, ...],
    algorithm_overrides: Path | None,
    defo_area_db: Path | None,
    event_db: Path | None,
    n_workers: int,
    threads_per_worker: int,
    block_shape: tuple[int, int],
    product_version: str,
    output_format: str,
    compression: bool,
    keep_relative: bool,
) -> None:
    r"""Create a calibration workflow configuration file.

    Generates a YAML configuration file with all parameters needed
    for displacement calibration, including input files, output paths,
    and processing settings.

    Examples
    --------
    Basic usage with required files:

        cal-disp config \
            --disp-file data/disp.h5 \
            --calibration-grid-latlon data/unr/grid_latlon_lookup_v0.2.txt \
            --calibration-grid-dir data/unr/ \
            --frame-id 8882 \
            --algorithm-params config/params.yaml \
            --los-file data/los.h5 \
            --dem-file data/dem.h5 \
            --output-dir outputs/ \
            --work-dir scratch/

    With optional corrections:

        cal-disp config \
            --disp-file data/disp.h5 \
            --calibration-grid-latlon data/unr/grid_latlon_lookup_v0.2.txt \
            --calibration-grid-dir data/unr/ \
            --frame-id 8882 \
            --algorithm-params config/params.yaml \
            --los-file data/los.h5 \
            --dem-file data/dem.h5 \
            --output-dir outputs/ \
            --work-dir scratch/ \
            --ref-tropo-files data/tropo_ref_1.h5 \
            --ref-tropo-files data/tropo_ref_2.h5 \
            --sec-tropo-files data/tropo_sec.h5 \
            --mask-file data/water_mask.tif

    """
    try:
        # Convert tuples to lists (or None if empty)
        ref_tropo_list = list(ref_tropo_files) if ref_tropo_files else None
        sec_tropo_list = list(sec_tropo_files) if sec_tropo_files else None
        iono_list = list(iono_files) if iono_files else None
        tiles_list = list(tiles_files) if tiles_files else None

        config_path = create_config(
            disp_file=disp_file,
            calibration_grid_latlon_file=calibration_grid_latlon,
            calibration_grid_ts_dir=calibration_grid_dir,
            frame_id=frame_id,
            algorithm_params_file=algorithm_params,
            los_file=los_file,
            dem_file=dem_file,
            output_dir=output_dir,
            work_dir=work_dir,
            config_name=config_file.name,
            mask_file=mask_file,
            ref_tropo_files=ref_tropo_list,
            sec_tropo_files=sec_tropo_list,
            iono_files=iono_list,
            tiles_files=tiles_list,
            algorithm_overrides_json=algorithm_overrides,
            defo_area_db_json=defo_area_db,
            event_db_json=event_db,
            n_workers=n_workers,
            threads_per_worker=threads_per_worker,
            block_shape=block_shape,
            product_version=product_version,
            output_format=output_format,
            compression=compression,
            keep_relative=keep_relative,
        )
        click.echo(f"Configuration created: {config_path}")

    except (FileNotFoundError, ValueError) as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()
