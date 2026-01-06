from __future__ import annotations

import glob
import json
from pathlib import Path
from typing import Annotated, Any

from pydantic import BeforeValidator


def _read_file_list_or_glob(_cls, value):
    """Convert input to list of file paths.

    Auto-detects:
    - Glob patterns: "*.nc"
    - File lists: text files containing paths (one per line)
    - Single files: anything else
    """
    if value is None:
        return []

    if isinstance(value, dict):
        raise TypeError(f"Expected string, Path, or list, but got dict: {value}")

    # Handle lists
    if isinstance(value, (list, tuple)):
        if len(value) == 1 and glob.has_magic(str(value[0])):
            return [Path(f) for f in glob.glob(str(value[0]))]
        return [Path(f) for f in value]

    # Handle string or Path
    v_path = Path(value)

    # Glob pattern
    if glob.has_magic(str(value)):
        return [Path(f) for f in glob.glob(str(value))]

    # Try to parse as file list
    if v_path.exists() and v_path.is_file():
        try:
            lines = [
                line.strip() for line in v_path.read_text().splitlines() if line.strip()
            ]

            # Empty file or binary → treat as single file
            if not lines:
                return [v_path]

            # Check if it looks like a file list by seeing if referenced files exist
            parent = v_path.parent
            resolved_paths = [
                parent / Path(f) if not Path(f).is_absolute() else Path(f)
                for f in lines
            ]

            # If at least one referenced path exists, treat as file list
            if any(p.exists() for p in resolved_paths):
                return resolved_paths

        except (UnicodeDecodeError, OSError):
            # Can't read as text → treat as single file
            pass

    # Single file (whether it exists or not)
    return [v_path]


def _to_path_required(v: str | Path | None) -> Path:
    """Convert to Path, rejecting None/empty."""
    if v is None:
        raise ValueError("Path cannot be None")
    if isinstance(v, str):
        if not v.strip():
            raise ValueError("Path cannot be an empty string")
        return Path(v)
    return v


def _to_path_optional(v: str | Path | None) -> Path | None:
    """Convert to Path, allowing None."""
    if v is None:
        return None
    if isinstance(v, str):
        if not v.strip():
            return None
        return Path(v)
    return v


def _to_path_or_cwd(v: str | Path | None) -> Path:
    """Convert to Path, using cwd for None/empty."""
    if v is None or v == "":
        return Path()
    return Path(v) if isinstance(v, str) else v


def _to_existing_file(v: str | Path | None) -> Path:
    """Convert to Path and verify it exists."""
    p = _to_path_required(v)
    if not p.exists():
        raise ValueError(f"File does not exist: {p}")
    return p


# Type aliases with validation
RequiredPath = Annotated[Path, BeforeValidator(_to_path_required)]
OptionalPath = Annotated[Path | None, BeforeValidator(_to_path_optional)]
DirectoryPath = Annotated[Path, BeforeValidator(_to_path_or_cwd)]
ExistingFilePath = Annotated[Path, BeforeValidator(_to_existing_file)]


def convert_paths_to_strings(obj: Any) -> Any:
    """Recursively convert Path objects to strings."""
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, dict):
        return {k: convert_paths_to_strings(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [convert_paths_to_strings(item) for item in obj]
    return obj


def format_summary_section(
    title: str, items: dict[str, Any], max_width: int = 70
) -> list[str]:
    """Format a section for summary output."""
    lines = [title, "=" * max_width, ""]
    for key, value in items.items():
        lines.append(f"  {key}: {value}")
    return lines


def _parse_algorithm_overrides(
    overrides_file: Path | str | None, frame_id: int | str
) -> dict[str, Any]:
    """Find frame-specific parameters to override for algorithm_parameters."""
    if overrides_file is None:
        return {}

    with open(overrides_file) as f:
        overrides = json.load(f)
        # Handle both {"data": {frame_id: ...}} and {frame_id: ...} formats
        if "data" in overrides:
            return overrides["data"].get(str(frame_id), {})
        return overrides.get(str(frame_id), {})
