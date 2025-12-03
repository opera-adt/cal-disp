from __future__ import annotations

from pathlib import Path
from typing import Any, ClassVar, Dict, Optional

from pydantic import ConfigDict, Field, field_validator

from ._yaml import YamlModel
from .workflow import (
    CalibrationWorkflow,
    DynamicAncillaryFileGroup,
    InputFileGroup,
    StaticAncillaryFileGroup,
    WorkerSettings,
)


class PrimaryExecutable(YamlModel):
    """Group describing the primary executable.

    Attributes
    ----------
    product_type : str
        Product type identifier for the PGE.

    """

    product_type: str = Field(
        default="CAL_DISP",
        description="Product type of the PGE.",
    )

    model_config = ConfigDict(extra="forbid")


class OutputOptions(YamlModel):
    """Output configuration options.

    Attributes
    ----------
    product_version : str
        Version of the product in <major>.<minor> format.
    output_format : str
        Format for output files (e.g., 'netcdf', 'hdf5').
    compression : bool
        Whether to compress output files.

    """

    product_version: str = Field(
        default="1.0",
        description="Version of the product, in <major>.<minor> format.",
    )

    output_format: str = Field(
        default="netcdf",
        description="Output file format.",
    )

    compression: bool = Field(
        default=True,
        description="Whether to compress output files.",
    )

    model_config = ConfigDict(extra="forbid")


class ProductPathGroup(YamlModel):
    """Group describing the product paths.

    Attributes
    ----------
    product_path : Path
        Directory where PGE will place results.
    scratch_path : Path
        Path to the scratch directory for intermediate files.
    output_path : Path
        Path to the SAS output directory.

    """

    product_path: Path = Field(
        default=Path(),
        description="Directory where PGE will place results.",
    )

    scratch_path: Path = Field(
        default=Path("./scratch"),
        description="Path to the scratch directory for intermediate files.",
    )

    output_path: Path = Field(
        default=Path("./output"),
        description="Path to the SAS output directory.",
        alias="sas_output_path",
    )

    @field_validator("product_path", "scratch_path", "output_path", mode="before")
    @classmethod
    def _validate_paths(cls, v: Any) -> Path:
        """Validate and convert directory paths."""
        if v is None or v == "":
            return Path()
        if isinstance(v, str):
            return Path(v)
        return v

    model_config = ConfigDict(extra="forbid", populate_by_name=True)


class RunConfig(YamlModel):
    """A PGE (Product Generation Executive) run configuration.

    This class represents the top-level configuration for running the
    calibration workflow as a PGE. It includes input files, output options,
    paths, and worker settings.

    Attributes
    ----------
    input_file_group : InputFileGroup
        Configuration for input files.
    dynamic_ancillary_file_group : Optional[DynamicAncillaryFileGroup]
        Dynamic ancillary files configuration.
    static_ancillary_file_group : Optional[StaticAncillaryFileGroup]
        Static ancillary files configuration.
    output_options : OutputOptions
        Output configuration options.
    primary_executable : PrimaryExecutable
        Primary executable configuration.
    product_path_group : ProductPathGroup
        Product path configuration.
    worker_settings : WorkerSettings
        Dask worker and parallelism configuration.
    log_file : Optional[Path]
        Path to the output log file.

    """

    # Used for the top-level key in YAML
    name: ClassVar[str] = "cal_disp_workflow"

    # Input configuration
    input_file_group: InputFileGroup = Field(
        ...,
        description="Configuration for required input files.",
    )

    dynamic_ancillary_file_group: Optional[DynamicAncillaryFileGroup] = Field(
        default=None,
        description="Dynamic ancillary files configuration.",
    )

    static_ancillary_file_group: Optional[StaticAncillaryFileGroup] = Field(
        default=None,
        description="Static ancillary files configuration.",
    )

    # Output and execution configuration
    output_options: OutputOptions = Field(
        default_factory=OutputOptions,
        description="Output configuration options.",
    )

    primary_executable: PrimaryExecutable = Field(
        default_factory=PrimaryExecutable,
        description="Primary executable configuration.",
    )

    product_path_group: ProductPathGroup = Field(
        default_factory=ProductPathGroup,
        description="Product path configuration.",
    )

    # Worker settings
    worker_settings: WorkerSettings = Field(
        default_factory=WorkerSettings,
        description="Dask worker and parallelism configuration.",
    )

    # Logging
    log_file: Optional[Path] = Field(
        default=None,
        description=(
            "Path to the output log file in addition to logging to stderr. "
            "If None, will be set based on output path."
        ),
    )

    def to_workflow(self) -> CalibrationWorkflow:
        """Convert PGE RunConfig to a CalibrationWorkflow object.

        This method translates the PGE-style configuration into the format
        expected by CalibrationWorkflow.

        Returns
        -------
        CalibrationWorkflow
            Converted workflow configuration.

        Examples
        --------
        >>> run_config = RunConfig.from_yaml_file("pge_config.yaml")
        >>> workflow = run_config.to_workflow()
        >>> workflow.create_directories()
        >>> workflow.run()

        """
        # Set up directories
        scratch_directory = self.product_path_group.scratch_path
        output_directory = self.product_path_group.output_path

        # Set up log file
        log_file = self.log_file
        if log_file is None:
            log_file = output_directory / "cal_disp_workflow.log"

        # Create the workflow
        workflow = CalibrationWorkflow(
            # Input files
            input_options=self.input_file_group,
            # Ancillary files
            dynamic_ancillary_file_options=self.dynamic_ancillary_file_group,
            static_ancillary_file_options=self.static_ancillary_file_group,
            # Directories
            work_directory=scratch_directory,
            output_directory=output_directory,
            # Settings
            worker_settings=self.worker_settings,
            log_file=log_file,
            # Don't resolve paths yet - let the workflow handle it
            keep_paths_relative=False,
        )

        return workflow

    def create_directories(self, exist_ok: bool = True) -> None:
        """Create all necessary directories for the PGE run.

        Parameters
        ----------
        exist_ok : bool, default=True
            If True, don't raise error if directories already exist.

        """
        self.product_path_group.product_path.mkdir(parents=True, exist_ok=exist_ok)
        self.product_path_group.scratch_path.mkdir(parents=True, exist_ok=exist_ok)
        self.product_path_group.output_path.mkdir(parents=True, exist_ok=exist_ok)

        # Create tmp directory in scratch
        tmp_dir = self.product_path_group.scratch_path / "tmp"
        tmp_dir.mkdir(parents=True, exist_ok=exist_ok)

    def validate_ready_to_run(self) -> Dict[str, Any]:
        """Check if run configuration is ready for execution.

        Returns
        -------
        dict
            Validation results with 'ready', 'errors', and 'warnings' keys.

        """
        errors = []
        warnings = []

        # Check if input files are configured
        if not hasattr(self.input_file_group, "disp_file"):
            errors.append("input_file_group must have disp_file configured")

        if not hasattr(self.input_file_group, "calibration_reference_grid"):
            errors.append(
                "input_file_group must have calibration_reference_grid configured"
            )

        # Check dynamic ancillary files if required
        if self.dynamic_ancillary_file_group:
            if not hasattr(
                self.dynamic_ancillary_file_group, "algorithm_parameters_file"
            ):
                warnings.append(
                    "dynamic_ancillary_file_group missing algorithm_parameters_file"
                )
            if not hasattr(self.dynamic_ancillary_file_group, "geometry_file"):
                warnings.append("dynamic_ancillary_file_group missing geometry_file")

        return {
            "ready": len(errors) == 0,
            "errors": errors,
            "warnings": warnings,
        }

    def summary(self) -> str:
        """Generate a human-readable summary of the run configuration.

        Returns
        -------
        str
            Multi-line summary string.

        """
        lines = [
            "PGE Run Configuration",
            "=" * 70,
            "",
            "Product Information:",
            f"  Product type:     {self.primary_executable.product_type}",
            f"  Product version:  {self.output_options.product_version}",
            f"  Output format:    {self.output_options.output_format}",
            "",
            "Paths:",
            f"  Product path:     {self.product_path_group.product_path}",
            f"  Scratch path:     {self.product_path_group.scratch_path}",
            f"  Output path:      {self.product_path_group.output_path}",
            f"  Log file:         {self.log_file or 'Default'}",
            "",
            "Input Files:",
            f"  DISP file:        {self.input_file_group.disp_file}",
            f"  Calibration grid: {self.input_file_group.calibration_reference_grid}",
            "",
            "Worker Settings:",
            f"  Workers:          {self.worker_settings.n_workers}",
            f"  Threads/worker:   {self.worker_settings.threads_per_worker}",
            f"  Total threads:    {self.worker_settings.total_threads}",
        ]

        # Dynamic ancillary files
        if self.dynamic_ancillary_file_group:
            lines.extend(
                [
                    "",
                    "Dynamic Ancillary Files:",
                ]
            )
            dynamic_files = self.dynamic_ancillary_file_group.get_all_files()
            for name, path in list(dynamic_files.items())[:5]:
                lines.append(f"  {name}: {path}")
            if len(dynamic_files) > 5:
                lines.append(f"  ... and {len(dynamic_files) - 5} more")

        # Static ancillary files
        if self.static_ancillary_file_group:
            static_files = self.static_ancillary_file_group.get_all_files()
            if static_files:
                lines.extend(
                    [
                        "",
                        "Static Ancillary Files:",
                    ]
                )
                for name, path in static_files.items():
                    lines.append(f"  {name}: {path}")

        # Validation status
        status = self.validate_ready_to_run()
        lines.append("")
        if not status["ready"]:
            lines.append("⚠️  Status: NOT READY")
            lines.append("Errors:")
            for error in status["errors"]:
                lines.append(f"  - {error}")
        else:
            lines.append("✓ Status: READY")

        if status["warnings"]:
            lines.append("")
            lines.append("Warnings:")
            for warning in status["warnings"]:
                lines.append(f"  ⚠️  {warning}")

        lines.append("=" * 70)

        return "\n".join(lines)

    @classmethod
    def create_example(cls) -> "RunConfig":
        """Create an example PGE run configuration.

        Returns
        -------
        RunConfig
            Example configuration with placeholder values.

        """
        return cls(
            input_file_group=InputFileGroup(
                disp_file=Path("input/disp.h5"),
                calibration_reference_grid=Path("input/cal_grid.parquet"),
                frame_id=1,
            ),
            dynamic_ancillary_file_group=DynamicAncillaryFileGroup(
                algorithm_parameters_file=Path("config/algorithm_params.yaml"),
                geometry_file=Path("input/geometry.h5"),
            ),
            product_path_group=ProductPathGroup(
                scratch_path=Path("./scratch"), output_path=Path("./output")
            ),
            worker_settings=WorkerSettings.create_standard(),
        )

    @classmethod
    def from_yaml_file(cls, yaml_path: Path) -> "RunConfig":
        """Load run configuration from YAML file.

        Parameters
        ----------
        yaml_path : Path
            Path to YAML configuration file.

        Returns
        -------
        RunConfig
            Loaded and validated run configuration.

        """
        import yaml  # type: ignore[import-untyped]

        with open(yaml_path, "r") as f:
            data = yaml.safe_load(f)

        # Handle nested structure if present
        if cls.name in data:
            data = data[cls.name]

        return cls.model_validate(data)

    def to_yaml_file(self, yaml_path: Path, include_name_key: bool = True) -> None:
        """Save run configuration to YAML file.

        Parameters
        ----------
        yaml_path : Path
            Path where YAML file should be saved.
        include_name_key : bool, default=True
            If True, wrap the config in a dictionary with the name as key.

        """
        import yaml

        # Convert to dict
        data = self.model_dump(mode="python", by_alias=True)

        # Convert Path objects to strings
        def convert_paths(obj):
            if isinstance(obj, Path):
                return str(obj)
            elif isinstance(obj, dict):
                return {k: convert_paths(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_paths(item) for item in obj]
            return obj

        data = convert_paths(data)

        # Wrap in name key if requested
        if include_name_key:
            data = {self.name: data}

        yaml_path.parent.mkdir(parents=True, exist_ok=True)
        with open(yaml_path, "w") as f:
            yaml.dump(data, f, default_flow_style=False, sort_keys=False)

    model_config = ConfigDict(
        extra="forbid",
        validate_assignment=True,
        arbitrary_types_allowed=True,
    )
