from __future__ import annotations

import logging
from pathlib import Path
from typing import ClassVar, List, Optional

from pydantic import ConfigDict, Field, field_validator

from ._utils import DirectoryPath
from ._yaml import STRICT_CONFIG_WITH_ALIASES, ValidationResult, YamlModel
from .workflow import (
    CalibrationWorkflow,
    DynamicAncillaryFileGroup,
    InputFileGroup,
    StaticAncillaryFileGroup,
    WorkerSettings,
)

logger = logging.getLogger(__name__)

#: The only product type and output format this release writes.
PRODUCT_TYPE = "DISP_CAL"
OUTPUT_FORMAT = "netcdf"


class PrimaryExecutable(YamlModel):
    """Group describing the primary executable.

    Attributes
    ----------
    product_type : str
        Product type identifier for the PGE. Informational: the SAS always
        writes DISP-CAL products, and `RunConfig.to_workflow` warns when
        another value is given.

    """

    product_type: str = Field(
        default=PRODUCT_TYPE,
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
        Format for output files. Only 'netcdf' is supported in this release;
        any other value raises.
    compression : bool
        Whether to compress output files.
    processing_facility : str
        Product processing facility written to the product identification.
    product_data_access : str
        URL (or DOI) where the DISP-CAL products can be retrieved.
    static_layers_data_access : Optional[str]
        URL of the frame's DISP static layers product; None takes the value
        from the input DISP product.
    source_data_access : Optional[str]
        URL (or DOI) of the input DISP products; None takes the value from the
        input DISP product.

    """

    product_version: str = Field(
        default="1.0",
        description="Version of the product, in <major>.<minor> format.",
    )

    output_format: str = Field(
        default=OUTPUT_FORMAT,
        description="Output file format. Only 'netcdf' is supported.",
    )

    compression: bool = Field(
        default=True,
        description=(
            "Whether to compress output files (gzip level 4 + shuffle, (256, 256)"
            " chunks, as the DISP-S1 input)."
        ),
    )

    processing_facility: str = Field(
        default="NASA Jet Propulsion Laboratory on AWS",
        description="Product processing facility written to /identification.",
    )

    product_data_access: str = Field(
        default=(
            "https://search.asf.alaska.edu/#/?dataset=OPERA-S1&productTypes=DISP-S1-CAL"
        ),
        description="URL (or DOI) where this product can be retrieved.",
    )

    static_layers_data_access: Optional[str] = Field(
        default=None,
        description=(
            "URL of the DISP static layers product of the frame. If null, taken"
            " from the input DISP product's identification group."
        ),
    )

    source_data_access: Optional[str] = Field(
        default=None,
        description=(
            "URL (or DOI) where the input DISP products can be retrieved. If null,"
            " taken from the input DISP product's identification group."
        ),
    )

    @field_validator("output_format")
    @classmethod
    def _only_netcdf(cls, v: str) -> str:
        """Fail for formats the writer cannot produce."""
        if v != OUTPUT_FORMAT:
            raise ValueError(
                f"output_format={v!r} is not supported in this release; the product"
                f" is always written as {OUTPUT_FORMAT!r}"
            )
        return v

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

    product_path: DirectoryPath = Field(
        default=Path(),
        description="Directory where PGE will place results.",
    )

    scratch_path: DirectoryPath = Field(
        default=Path("./scratch"),
        description="Path to the scratch directory for intermediate files.",
    )

    output_path: DirectoryPath = Field(
        default=Path("./output"),
        description="Path to the SAS output directory.",
        alias="sas_output_path",
    )

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
    _yaml_root_key: ClassVar[str] = "cal_disp_workflow"

    # Input configuration
    input_file_group: InputFileGroup = Field(
        ...,
        description="Configuration for required input files.",
    )

    dynamic_ancillary_group: Optional[DynamicAncillaryFileGroup] = Field(
        default=None,
        description="Dynamic ancillary files configuration.",
    )

    static_ancillary_group: Optional[StaticAncillaryFileGroup] = Field(
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
        expected by CalibrationWorkflow. PGE fields the SAS does not use are
        reported with a WARNING instead of being dropped silently:
        ``product_path_group.product_path`` (products are written to
        ``sas_output_path``) and ``primary_executable.product_type``.

        Returns
        -------
        CalibrationWorkflow
            Converted workflow configuration.

        Examples
        --------
        >>> run_config = RunConfig.from_yaml_file("pge_config.yaml")
        >>> workflow = run_config.to_workflow()
        >>> workflow.create_directories()

        """
        paths = self.product_path_group
        if paths.product_path.resolve() != paths.output_path.resolve():
            logger.warning(
                "product_path_group.product_path (%s) is not used by the SAS:"
                " products are written to sas_output_path (%s)",
                paths.product_path,
                paths.output_path,
            )
        if self.primary_executable.product_type != PRODUCT_TYPE:
            logger.warning(
                "primary_executable.product_type=%r is not used: this SAS always"
                " writes %s products",
                self.primary_executable.product_type,
                PRODUCT_TYPE,
            )

        # Set up log file
        log_file = self.log_file
        if log_file is None:
            log_file = self.product_path_group.output_path / "cal_disp.log"

        # Create the workflow
        workflow = CalibrationWorkflow(
            # Input files
            input_options=self.input_file_group,
            # Ancillary files
            dynamic_ancillary_options=self.dynamic_ancillary_group,
            static_ancillary_options=self.static_ancillary_group,
            # Directories
            work_directory=self.product_path_group.scratch_path,
            output_directory=self.product_path_group.output_path,
            # Settings
            worker_settings=self.worker_settings,
            log_file=log_file,
            product_version=self.output_options.product_version,
            compression=self.output_options.compression,
            processing_facility=self.output_options.processing_facility,
            product_data_access=self.output_options.product_data_access,
            static_layers_data_access=self.output_options.static_layers_data_access,
            source_data_access=self.output_options.source_data_access,
            # Resolve paths to absolute
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

    def validate_ready_to_run(self) -> ValidationResult:
        """Check if run configuration is ready for execution."""
        errors: List[str] = []
        warnings: List[str] = []

        if self.input_file_group.disp_file is None:
            errors.append("disp_file must be provided")

        if self.input_file_group.unr_grid_latlon_file is None:
            errors.append("unr_grid_latlon_file must be provided")

        if self.input_file_group.unr_timeseries_dir is None:
            errors.append("unr_timeseries_dir must be provided")

        # Check for missing files only if required options are set
        if self.input_file_group is not None:
            missing = self.get_missing_files()
            if missing:
                warnings.append(f"Missing files: {', '.join(missing)}")

        # Check dynamic ancillaries
        if self.dynamic_ancillary_group is None:
            errors.append("dynamic_ancillary_options must be provided")

        # Check dynamic ancillary files if provided
        if self.dynamic_ancillary_group:
            if self.dynamic_ancillary_group.algorithm_parameters_file is None:
                warnings.append("Missing algorithm_parameters_file")
            if self.dynamic_ancillary_group.dem_file is None:
                warnings.append("Missing dem")
            if self.dynamic_ancillary_group.los_file is None:
                warnings.append("Missing los_file")

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
            f"  Frame ID:         {self.input_file_group.frame_id}",
            f"  UNR lookup:       {self.input_file_group.unr_grid_latlon_file}",
            f"  UNR grid dir:     {self.input_file_group.unr_timeseries_dir}",
            f"  UNR version:      {self.input_file_group.unr_grid_version}",
            f"  UNR type:         {self.input_file_group.unr_grid_type}",
            "",
            "Worker Settings:",
            f"  Workers:          {self.worker_settings.n_workers}",
            f"  Threads/worker:   {self.worker_settings.threads_per_worker}",
            f"  Total threads:    {self.worker_settings.total_threads}",
        ]

        # Dynamic ancillary files
        if self.dynamic_ancillary_group:
            lines.extend(
                [
                    "",
                    "Dynamic Ancillary Files:",
                ]
            )
            dynamic_files = self.dynamic_ancillary_group.get_all_files()
            for name, path in list(dynamic_files.items())[:5]:
                lines.append(f"  {name}: {path}")
            if len(dynamic_files) > 5:
                lines.append(f"  ... and {len(dynamic_files) - 5} more")

        # Static ancillary files
        if self.static_ancillary_group:
            static_files = self.static_ancillary_group.get_all_files()
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
            lines.append("Status: NOT READY")
            lines.append("Errors:")
            for error in status["errors"]:
                lines.append(f"  - {error}")
        else:
            lines.append("Status: READY")

        if status["warnings"]:
            lines.append("")
            lines.append("Warnings:")
            for warning in status["warnings"]:
                lines.append(f"  WARNING: {warning}")

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
        return cls.model_construct(
            input_file_group=InputFileGroup(
                disp_file=Path("input/disp_frame_8882.h5"),
                unr_grid_latlon_file=Path("input/unr/grid_latlon_lookup_v0.2.txt"),
                unr_timeseries_dir=Path("input/unr/"),
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
                skip_file_checks=True,
            ),
            dynamic_ancillary_group=DynamicAncillaryFileGroup(
                algorithm_parameters_file=Path("config/algorithm_params.yaml"),
                static_dem_file=Path("input/dem.tif"),
                static_los_file=Path("input/los.tif"),
            ),
            product_path_group=ProductPathGroup(
                scratch_path=Path("./scratch"), product_path=Path("./output")
            ),
            worker_settings=WorkerSettings.create_standard(),
        )

    @classmethod
    def from_yaml_file(cls, yaml_path: Path) -> "RunConfig":
        """Load run configuration from YAML file.

        Handles optional name wrapper key.

        Parameters
        ----------
        yaml_path : Path
            Path to YAML configuration file.

        Returns
        -------
        RunConfig
            Loaded and validated run configuration.

        """
        data = cls._load_yaml_data(yaml_path)

        # Get the root key safely (e.g., "cal_disp_workflow")
        root_key = getattr(cls, "_yaml_root_key", None)

        # Handle optional wrapper key if it exists in the data
        if root_key and isinstance(data, dict) and root_key in data:
            data = data[root_key]

        return cls.model_validate(data)

    @staticmethod
    def _load_yaml_data(yaml_path: Path) -> dict:
        """Load YAML data from file.

        Parameters
        ----------
        yaml_path : Path
            Path to YAML file.

        Returns
        -------
        dict
            Loaded YAML data.

        """
        from ruamel.yaml import YAML

        y = YAML(typ="safe")
        with open(yaml_path) as f:
            return y.load(f)

    model_config = STRICT_CONFIG_WITH_ALIASES
