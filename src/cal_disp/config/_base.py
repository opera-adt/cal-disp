from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from pydantic import ConfigDict, Field, field_validator

from ._utils import OptionalPath, RequiredPath, _read_file_list_or_glob
from ._yaml import STRICT_CONFIG_WITH_ALIASES, YamlModel


class InputFileGroup(YamlModel):
    """Input file group for the SAS.

    Attributes
    ----------
    disp_file : Path
        Path to DISP file.
    calibration_reference_grid : Path
        Path to UNR calibration reference file (parquet format).

    """

    disp_file: RequiredPath = Field(
        ...,
        description="Path to DISP file.",
    )

    calibration_reference_grid: RequiredPath = Field(
        ...,
        description="Path to UNR calibration reference file [parquet].",
    )

    frame_id: int = Field(
        ...,
        description="Frame ID of the DISP frame.",
    )

    model_config = ConfigDict(
        extra="forbid",
    )


class DynamicAncillaryFileGroup(YamlModel):
    """Dynamic ancillary files for the SAS.

    Attributes
    ----------
    algorithm_parameters_file : Path
        Path to file containing SAS algorithm parameters.
    geometry_file : Path
        Path to the DISP static_layer file with line-of-sight unit vectors.
    mask_file : Optional[Path]
        Optional byte mask file to ignore low correlation/bad data.
    troposphere_files : List[Path]
        Paths to TROPO files for atmospheric correction.
    ionosphere_files : List[Path]
        Paths to IONO files for ionospheric correction.
    tiles_files : Optional[List[Path]]
        Paths to calibration tile bounds files.

    """

    algorithm_parameters_file: RequiredPath = Field(
        ...,
        description="Path to file containing SAS algorithm parameters.",
    )

    geometry_file: RequiredPath = Field(
        ...,
        alias="static_layers_file",
        description=(
            "Path to the DISP static_layer file (1 per frame) with line-of-sight"
            " unit vectors."
        ),
    )

    mask_file: OptionalPath = Field(
        default=None,
        description=(
            "Optional Byte mask file used to ignore low correlation/bad data (e.g water"
            " mask). Convention is 0 for no data/invalid, and 1 for good data. Dtype"
            " must be uint8."
        ),
    )

    troposphere_files: List[Path] = Field(
        default_factory=list,
        alias="tropo_files",
        description=(
            "Paths to the TROPO files (1 for reference and 1 for secondary date)."
            " If none provided, tropospheric correction in calibration is skipped."
        ),
    )

    ionosphere_files: List[Path] = Field(
        default_factory=list,
        alias="iono_files",
        description=(
            "Paths to the IONO files (1 for reference and 1 for secondary date)."
            " If none provided, ionospheric correction in calibration is skipped."
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
        "troposphere_files", "ionosphere_files", "tiles_files", mode="before"
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
