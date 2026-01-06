from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from pydantic import ConfigDict, Field, field_validator, model_validator

from ._utils import OptionalPath, RequiredPath, _read_file_list_or_glob
from ._yaml import STRICT_CONFIG_WITH_ALIASES, YamlModel


class InputFileGroup(YamlModel):
    """Input file group for the SAS.

    Attributes
    ----------
    disp_file : Path
        Path to DISP file.
    calibration_reference_latlon_file : Path
        Path to UNR grid lookup table (e.g., grid_latlon_lookup_v0.2.txt).
    calibration_reference_grid_dir : Path
        Directory containing UNR .tenv8 timeseries files.
    frame_id : int
        Frame ID of the DISP frame.

    """

    disp_file: RequiredPath = Field(
        ...,
        description="Path to DISP file.",
    )

    calibration_reference_latlon_file: RequiredPath = Field(
        ...,
        description="Path to UNR grid lookup table (grid_latlon_lookup_v0.2.txt).",
    )

    calibration_reference_grid_dir: RequiredPath = Field(
        ...,
        description="Directory containing UNR .tenv8 timeseries files.",
    )

    frame_id: int = Field(
        ...,
        description="Frame ID of the DISP frame.",
    )

    skip_file_checks: bool = False

    model_config = ConfigDict(
        extra="forbid",
    )

    @field_validator("disp_file")
    @classmethod
    def validate_disp_file(cls, v: Path) -> Path:
        """Validate DISP file has .nc extension."""
        if v.suffix not in [".nc", ".h5"]:
            msg = f"DISP file must be .nc or .h5, got {v.suffix}"
            raise ValueError(msg)
        return v

    @field_validator("calibration_reference_latlon_file")
    @classmethod
    def validate_latlon_file(cls, v: Path) -> Path:
        """Validate latlon file is a lookup table."""
        if not v.name.startswith("grid_latlon_lookup"):
            msg = f"Expected grid_latlon_lookup file, got {v.name}"
            raise ValueError(msg)
        if not v.name.endswith(".txt"):
            msg = f"Lookup table must be .txt, got {v.suffix}"
            raise ValueError(msg)
        return v

    @field_validator("frame_id")
    @classmethod
    def validate_frame_id(cls, v: int) -> int:
        """Validate frame ID is reasonable."""
        if not 1 <= v <= 99999:
            msg = f"Frame ID must be between 1 and 99999, got {v}"
            raise ValueError(msg)
        return v

    @model_validator(mode="after")
    def validate_grid_dir(self) -> "InputFileGroup":
        """Validate grid directory contains .tenv8 files."""
        if self.skip_file_checks:
            return self

        v = self.calibration_reference_grid_dir
        if not v.exists():
            raise ValueError(f"Grid directory does not exist: {v}")
        if not v.is_dir():
            raise ValueError(f"Grid directory is not a directory: {v}")
        if not any(v.glob("*.tenv8")):
            raise ValueError(f"No .tenv8 files found in {v}")

        return self


class DynamicAncillaryFileGroup(YamlModel):
    """Dynamic ancillary files for the SAS.

    Attributes
    ----------
    algorithm_parameters_file : Path
        Path to file containing SAS algorithm parameters.
    los_file : Path
        Path to the DISP static LOS layer file (line-of-sight unit vectors).
        Alias: static_los_file
    dem_file : Path
        Path to the DISP static DEM layer file (digital elevation model).
        Alias: static_dem_file
    mask_file : Path or None, optional
        Optional byte mask file to ignore low correlation/bad data (e.g., water mask).
        Convention: 0 = invalid/no data, 1 = good data. Dtype must be uint8.
        Default is None.
    reference_tropo_files : list[Path] or None, optional
        Paths to TROPO files for the reference (primary) date.
        If not provided, tropospheric correction for reference is skipped.
        Alias: ref_tropo_files. Default is None.
    secondary_tropo_files : list[Path] or None, optional
        Paths to TROPO files for the secondary date.
        If not provided, tropospheric correction for secondary is skipped.
        Alias: sec_tropo_files. Default is None.
    iono_files : list[Path] or None, optional
        Paths to ionospheric correction files.
        If not provided, ionospheric correction is skipped.
        Default is None.
    tiles_files : list[Path] or None, optional
        Paths to calibration tile bounds files (e.g., S1 burst bounds) covering
        the full frame. If not provided, per-tile calibration is skipped.
        Default is None.

    """

    algorithm_parameters_file: RequiredPath = Field(
        ...,
        description="Path to file containing SAS algorithm parameters.",
    )

    los_file: RequiredPath = Field(
        ...,
        alias="static_los_file",
        description=(
            "Path to the DISP static los layer file (1 per frame) with line-of-sight"
            " unit vectors."
        ),
    )

    dem_file: RequiredPath = Field(
        ...,
        alias="static_dem_file",
        description=(
            "Path to the DISP static dem layer file (1 per frame) with line-of-sight"
            " unit vectors."
        ),
    )
    # NOTE should I add also shadow_layover static file as input

    mask_file: OptionalPath = Field(
        default=None,
        description=(
            "Optional Byte mask file used to ignore low correlation/bad data (e.g water"
            " mask). Convention is 0 for no data/invalid, and 1 for good data. Dtype"
            " must be uint8."
        ),
    )

    reference_tropo_files: Optional[List[Path]] = Field(
        default=None,
        alias="ref_tropo_files",
        description=(
            "Path to the TROPO file for the reference date."
            " If not provided, tropospheric correction for reference is skipped."
        ),
    )

    secondary_tropo_files: Optional[List[Path]] = Field(
        default=None,
        alias="sec_tropo_files",
        description=(
            "Path to the TROPO file for the secondary date."
            " If not provided, tropospheric correction for secondary is skipped."
        ),
    )

    iono_files: Optional[List[Path]] = Field(
        default=None,
        alias="iono_files",
        description=(
            "Path to the IONO files"
            " If not provided, ionosphere correction for reference is skipped."
        ),
    )

    tiles_files: Optional[List[Path]] = Field(
        default=None,
        description=(
            "Paths to the calibration tile bounds files (e.g. S1 burst bounds) covering"
            " full frame. If none provided, calibration per tile is skipped."
        ),
    )

    @field_validator(
        "reference_tropo_files",
        "secondary_tropo_files",
        "iono_files",
        "tiles_files",
        mode="before",
    )
    @classmethod
    def _validate_file_lists(cls, v):
        """Validate and process file lists or glob patterns."""
        return _read_file_list_or_glob(cls, v)

    def get_all_files(self) -> Dict[str, Path | list[Path]]:
        return self.get_all_file_paths(flatten_lists=True)

    model_config = STRICT_CONFIG_WITH_ALIASES


class StaticAncillaryFileGroup(YamlModel):
    """Static ancillary files for the SAS.

    These files contain configuration and reference data that don't change
    between processing runs for a given frame.

    Attributes
    ----------
    algorithm_parameters_overrides_json : Optional[Path]
        JSON file with frame-specific algorithm parameter overrides.
    deformation_area_database_json : Optional[Path]
        GeoJSON file with deforming areas to exclude from calibration.
    event_database_json : Optional[Path]
        GeoJSON file with earthquake/volcanic activity events for each frame.

    """

    algorithm_parameters_overrides_json: OptionalPath = Field(
        default=None,
        description=(
            "JSON file containing frame-specific algorithm parameters to override the"
            " defaults passed in the `algorithm_parameters.yaml`."
        ),
    )

    deformation_area_database_json: OptionalPath = Field(
        default=None,
        alias="defo_area_db_json",
        description=(
            "GeoJSON file containing list of deforming areas to exclude from"
            " calibration (e.g. Central Valley subsidence)."
        ),
    )

    event_database_json: OptionalPath = Field(
        default=None,
        alias="event_db_json",
        description=(
            "GeoJSON file containing list of events (earthquakes, volcanic activity)"
            " for each frame."
        ),
    )

    def has_algorithm_overrides(self) -> bool:
        """Check if algorithm parameter overrides are provided."""
        return self.algorithm_parameters_overrides_json is not None

    def has_deformation_database(self) -> bool:
        """Check if deformation area database is provided."""
        return self.deformation_area_database_json is not None

    def has_event_database(self) -> bool:
        """Check if event database is provided."""
        return self.event_database_json is not None

    def get_all_files(self) -> Dict[str, Path | list[Path]]:
        return self.get_all_file_paths(flatten_lists=True)

    model_config = STRICT_CONFIG_WITH_ALIASES
