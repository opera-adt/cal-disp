from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

from pydantic import Field, model_validator

from ._base import DynamicAncillaryFileGroup, InputFileGroup, StaticAncillaryFileGroup
from ._utils import DirectoryPath
from ._workers import WorkerSettings
from ._yaml import STRICT_CONFIG_WITH_ALIASES, ValidationResult, YamlModel


class CalibrationWorkflow(YamlModel):
    """Calibration workflow configuration.

    This class manages the complete configuration for the displacement calibration
    workflow, including input files, output directories, worker settings, and
    logging configuration.

    Attributes
    ----------
    input_options : Optional[InputFileGroup]
        Configuration for required input files.
    work_directory : Path
        Directory for intermediate processing files.
    output_directory : Path
        Directory for final output files.
    keep_paths_relative : bool
        If False, resolve all paths to absolute paths.
    dynamic_ancillary_file_options : Optional[DynamicAncillaryFileGroup]
        Optional dynamic ancillary files for processing.
    static_ancillary_file_options : Optional[StaticAncillaryFileGroup]
        Optional static ancillary files for processing.
    worker_settings : WorkerSettings
        Dask worker and threading configuration.
    log_file : Optional[Path]
        Custom log file path.

    """

    # Input/output file configuration
    input_options: Optional[InputFileGroup] = Field(
        default=None,
        description=(
            "Configuration for required input files. Must be provided before running"
            " workflow."
        ),
    )

    work_directory: DirectoryPath = Field(
        default=Path(),
        description=(
            "Directory for intermediate processing files. Created if it doesn't exist."
        ),
    )

    output_directory: DirectoryPath = Field(
        default=Path(),
        description="Directory for final output files. Created if it doesn't exist.",
    )

    keep_paths_relative: bool = Field(
        default=False,
        description=(
            "If False, resolve all relative paths to absolute paths. "
            "If True, keep paths as provided."
        ),
    )

    # Optional ancillary file groups
    dynamic_ancillary_options: Optional[DynamicAncillaryFileGroup] = Field(
        default=None,
        description="Dynamic ancillary files (dem, los, masks, troposphere, etc.).",
    )

    static_ancillary_options: Optional[StaticAncillaryFileGroup] = Field(
        default=None,
        description="Static ancillary files (algorithm overrides, databases, etc.).",
    )

    # Worker and logging configuration
    worker_settings: WorkerSettings = Field(
        default_factory=WorkerSettings,
        description="Worker and parallelism configuration.",
    )

    log_file: Optional[Path] = Field(
        default=None,
        description=(
            "Path to output log file (in addition to logging to stderr). "
            "If None, defaults to 'cal_disp.log' in work_directory."
        ),
    )

    @model_validator(mode="after")
    def _resolve_paths(self) -> "CalibrationWorkflow":
        """Resolve all paths to absolute if keep_paths_relative is False.

        Uses object.__setattr__() to bypass validate_assignment and avoid recursion.
        """
        if not self.keep_paths_relative:
            object.__setattr__(self, "work_directory", self.work_directory.resolve())
            object.__setattr__(
                self, "output_directory", self.output_directory.resolve()
            )

            if self.log_file is not None:
                object.__setattr__(self, "log_file", self.log_file.resolve())

        return self

    @model_validator(mode="after")
    def _set_default_log_file(self) -> "CalibrationWorkflow":
        """Set default log file path if not provided."""
        if self.log_file is None:
            object.__setattr__(self, "log_file", self.work_directory / "cal_disp.log")

        return self

    def validate_ready_to_run(self) -> ValidationResult:
        """Check if workflow is ready to run."""
        errors = []
        warnings = []

        if self.input_options is None:
            errors.append("input_options must be provided")
        else:
            if self.input_options.disp_file is None:
                errors.append("disp_file must be provided in input_options")
            if self.input_options.calibration_reference_latlon_file is None:
                errors.append(
                    "calibration_reference_latlon_file must be provided in"
                    " input_options"
                )
            if self.input_options.calibration_reference_grid_dir is None:
                errors.append(
                    "calibration_reference_grid_dir must be provided in input_options"
                )

        # Check dynamic ancillaries
        if self.dynamic_ancillary_options is None:
            errors.append("dynamic_ancillary_options must be provided")

        # Check for missing files only if required options are set
        if self.input_options is not None:
            missing = self.get_missing_files()
            if missing:
                warnings.append(f"Missing files: {', '.join(missing)}")

        return {
            "ready": len(errors) == 0,
            "errors": errors,
            "warnings": warnings,
        }

    def validate_input_files_exist(self) -> Dict[str, Dict[str, Any]]:
        """Check if all input files exist."""
        results = {}

        if self.input_options:
            results.update(self.input_options.validate_files_exist())

        if self.dynamic_ancillary_options:
            results.update(self.dynamic_ancillary_options.validate_files_exist())

        if self.static_ancillary_options:
            results.update(self.static_ancillary_options.validate_files_exist())

        return results

    def get_missing_files(self) -> List[str]:
        """Get list of missing required input files."""
        status = self.validate_input_files_exist()
        return [name for name, info in status.items() if not info["exists"]]

    def create_directories(self, exist_ok: bool = True) -> None:
        """Create work and output directories if they don't exist.

        Parameters
        ----------
        exist_ok : bool, default=True
            If True, don't raise error if directories already exist.

        """
        self.work_directory.mkdir(parents=True, exist_ok=exist_ok)
        self.output_directory.mkdir(parents=True, exist_ok=exist_ok)

        # Also create parent directory for log file
        if self.log_file is not None:
            self.log_file.parent.mkdir(parents=True, exist_ok=exist_ok)

    def setup_logging(
        self, level: int = 20, format_string: Optional[str] = None  # logging.INFO
    ):
        """Set up logging configuration for the workflow.

        Parameters
        ----------
        level : int, default=20 (INFO)
            Logging level (DEBUG=10, INFO=20, WARNING=30, ERROR=40, CRITICAL=50).
        format_string : str, optional
            Custom format string for log messages.

        Returns
        -------
        logging.Logger
            Configured logger instance.

        """
        import logging

        if format_string is None:
            format_string = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

        logger = logging.getLogger("cal_disp")
        logger.setLevel(level)
        logger.handlers.clear()

        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(level)
        console_handler.setFormatter(logging.Formatter(format_string))
        logger.addHandler(console_handler)

        # File handler
        if self.log_file is not None:
            self.log_file.parent.mkdir(parents=True, exist_ok=True)
            file_handler = logging.FileHandler(self.log_file)
            file_handler.setLevel(level)
            file_handler.setFormatter(logging.Formatter(format_string))
            logger.addHandler(file_handler)

        return logger

    def summary(self) -> str:
        """Generate a human-readable summary of the workflow configuration.

        Returns
        -------
        str
            Multi-line summary string.

        """
        lines = [
            "Calibration Workflow Configuration",
            "=" * 70,
            "",
            "Directories:",
            f"  Work directory:   {self.work_directory}",
            f"  Output directory: {self.output_directory}",
            f"  Log file:         {self.log_file}",
            f"  Keep relative:    {self.keep_paths_relative}",
            "",
        ]

        # Input files
        if self.input_options:
            lines.extend(
                [
                    "Input Files:",
                    f"  DISP file:        {self.input_options.disp_file}",
                    f"  Frame ID:         {self.input_options.frame_id}",
                    (
                        "  UNR lookup:       "
                        f"{self.input_options.calibration_reference_latlon_file}"
                    ),
                    (
                        "  UNR grid dir:     "
                        f"{self.input_options.calibration_reference_grid_dir}"
                    ),
                    "",
                ]
            )
        else:
            lines.extend(
                [
                    "Input Files:",
                    "  Not configured!",
                    "",
                ]
            )

        # Worker settings
        lines.extend(
            [
                "Worker Settings:",
                f"  Workers:          {self.worker_settings.n_workers}",
                f"  Threads/worker:   {self.worker_settings.threads_per_worker}",
                f"  Total threads:    {self.worker_settings.total_threads}",
                f"  Block shape:      {self.worker_settings.block_shape}",
                "",
            ]
        )

        # Dynamic ancillary files
        if self.dynamic_ancillary_options:
            dynamic_files = self.dynamic_ancillary_options.get_all_files()
            if dynamic_files:
                lines.extend(["Dynamic Ancillary Files:"])
                for name, path in list(dynamic_files.items())[:5]:  # Show first 5
                    lines.append(f"  {name}: {path}")
                if len(dynamic_files) > 5:
                    lines.append(f"  ... and {len(dynamic_files) - 5} more")
                lines.append("")

        # Static ancillary files
        if self.static_ancillary_options:
            static_files = self.static_ancillary_options.get_all_files()
            if static_files:
                lines.extend(["Static Ancillary Files:"])
                for name, path in static_files.items():
                    lines.append(f"  {name}: {path}")
                lines.append("")

        # Check readiness
        status = self.validate_ready_to_run()
        if not status["ready"]:
            lines.extend(
                [
                    "Workflow Status: NOT READY",
                    "Errors:",
                ]
            )
            for error in status["errors"]:
                lines.append(f"  - {error}")
        else:
            lines.extend(
                [
                    "Workflow Status: READY",
                ]
            )

        if status["warnings"]:
            lines.extend(
                [
                    "",
                    "Warnings:",
                ]
            )
            for warning in status["warnings"]:
                lines.append(f"  WARNING: {warning}")

        lines.append("=" * 70)

        return "\n".join(lines)

    @classmethod
    def create_example(cls) -> "CalibrationWorkflow":
        """Create an example workflow configuration.

        Returns
        -------
        CalibrationWorkflow
            Example configuration with placeholder values.

        """
        return cls(
            work_directory=Path("./work"),
            output_directory=Path("./output"),
            input_options=InputFileGroup(
                disp_file=Path("input/disp_frame_8882.nc"),
                calibration_reference_latlon_file=Path(
                    "input/unr/grid_latlon_lookup_v0.2.txt"
                ),
                calibration_reference_grid_dir=Path("input/unr/"),
                frame_id=8882,
                skip_file_checks=True,
            ),
            dynamic_ancillary_options=DynamicAncillaryFileGroup(
                algorithm_parameters_file="algorithm.yaml",
                static_los_file="line_of_sight_enu.tif",
                static_dem_file="dem.tif",
            ),
            worker_settings=WorkerSettings.create_standard(),
            keep_paths_relative=True,
        )

    @classmethod
    def create_minimal(cls) -> "CalibrationWorkflow":
        """Create a minimal workflow configuration without input files.

        Returns
        -------
        CalibrationWorkflow
            Minimal configuration. Input files must be added before running.

        """
        return cls(
            work_directory=Path("./work"),
            output_directory=Path("./output"),
            keep_paths_relative=True,
        )

    model_config = STRICT_CONFIG_WITH_ALIASES
