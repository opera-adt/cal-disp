from __future__ import annotations

import glob
import json
from pathlib import Path
from typing import Annotated, Any, Dict, List, Optional, Union

from pydantic import BeforeValidator


def _read_file_list_or_glob(_cls, value):
    """Check if the input file list is a glob pattern or a text file.

    Parameters
    ----------
    _cls : type
        Pydantic model class
    value : str | Path | list[str] | list[Path] | None
        Input file list, glob pattern, or path to text file containing file list.

    Returns
    -------
    list[Path]
        List of Path objects representing files.

    Raises
    ------
    ValueError
        If the provided path does not exist or is not a file.
    TypeError
        If value is not a valid type (dict, etc.)

    """
    if value is None:
        return []

    # Reject invalid types
    if isinstance(value, dict):
        msg = f"Expected string, Path, or list, but got dict: {value}"
        raise TypeError(msg)

    # Check if they've passed a glob pattern
    if (
        isinstance(value, (list, tuple))
        and len(value) == 1
        and glob.has_magic(str(value[0]))
    ):
        value = glob.glob(str(value[0]))
    elif isinstance(value, (str, Path)):
        v_path = Path(value)
        # Check if it's a glob pattern
        if glob.has_magic(str(value)):
            value = glob.glob(str(value))
        # Check if it's a newline-delimited list of input files
        elif v_path.exists() and v_path.is_file():
            filenames = [Path(f) for f in v_path.read_text().splitlines() if f.strip()]
            parent = v_path.parent
            return [parent / f if not f.is_absolute() else f for f in filenames]
        else:
            # Don't raise error for non-existent files in list
            # Just treat it as a single file path
            return [v_path]

    return [Path(f) for f in value]


def validate_path_field(
    v: Union[str, Path, None], allow_none: bool = False, allow_empty: bool = False
) -> Optional[Path]:
    """Reusable path validator.

    Parameters
    ----------
    v : str | Path | None
        Value to validate.
    allow_none : bool, default=False
        Whether to allow None values.
    allow_empty : bool, default=False
        Whether to allow empty strings.

    Returns
    -------
    Path | None
        Validated path.

    Raises
    ------
    ValueError
        If validation fails.

    """
    if v is None:
        if allow_none:
            return None
        raise ValueError("Path cannot be None")

    if isinstance(v, str):
        if not v.strip():
            if allow_empty:
                return None
            raise ValueError("Path cannot be an empty string")
        return Path(v)

    return v


def _validate_directory_path(v: Union[str, Path, None]) -> Path:
    """Validate and convert to Path, allowing empty for current directory.

    Parameters
    ----------
    v : str | Path | None
        Value to validate.

    Returns
    -------
    Path
        Validated Path object (empty Path() if None or empty string).

    """
    if v is None or v == "":
        return Path()

    if isinstance(v, str):
        return Path(v)

    return v


# Create specific validators
def _to_path_required(v: Union[str, Path, None]) -> Path:
    if v is None:
        raise ValueError("Path cannot be None")

    if isinstance(v, str):
        if not v.strip():
            raise ValueError("Path cannot be an empty string")
        return Path(v)

    # v must be Path at this point
    return v


def _to_path_optional(v: Union[str, Path, None]) -> Optional[Path]:
    """Convert to Path, optional."""
    return validate_path_field(v, allow_none=True)


def _to_existing_file(v: Union[str, Path, None]) -> Path:
    p = _to_path_required(v)  # first do the basic checks
    if not p.exists():
        raise ValueError(f"File does not exist: {p}")
    return p


# Type aliases for cleaner code
RequiredPath = Annotated[Path, BeforeValidator(_to_path_required)]
OptionalPath = Annotated[Optional[Path], BeforeValidator(_to_path_optional)]
DirectoryPath = Annotated[Path, BeforeValidator(_validate_directory_path)]
ExistingFilePath = Annotated[Path, BeforeValidator(_to_existing_file)]


def convert_paths_to_strings(obj: Any) -> Any:
    """Recursively convert Path objects to strings in nested structures.

    Parameters
    ----------
    obj : Any
        Object potentially containing Path objects.

    Returns
    -------
    Any
        Same structure with Path objects converted to strings.

    """
    if isinstance(obj, Path):
        return str(obj)
    elif isinstance(obj, dict):
        return {k: convert_paths_to_strings(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_paths_to_strings(item) for item in obj]
    return obj


def format_summary_section(
    title: str, items: Dict[str, Any], max_width: int = 70
) -> List[str]:
    """Format a section for summary output."""
    lines = [title, "=" * max_width, ""]
    for key, value in items.items():
        lines.append(f"  {key}: {value}")
    return lines


# WORKFLOW UTILS
def _parse_algorithm_overrides(
    overrides_file: Path | str | None, frame_id: int | str
) -> dict[str, Any]:
    """Find the frame-specific parameters to override for algorithm_parameters."""
    if overrides_file is not None:
        with open(overrides_file) as f:
            overrides = json.load(f)
            if "data" in overrides:
                return overrides["data"].get(str(frame_id), {})
            else:
                return overrides.get(str(frame_id), {})
    return {}
