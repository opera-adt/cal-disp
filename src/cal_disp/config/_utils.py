from __future__ import annotations

import glob
import json
from pathlib import Path
from typing import Any


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
