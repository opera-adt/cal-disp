"""Tests for configuration utility functions."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cal_disp.config._utils import (
    _parse_algorithm_overrides,
    _read_file_list_or_glob,
    convert_paths_to_strings,
    format_summary_section,
)


class TestReadFileListOrGlob:
    """Tests for _read_file_list_or_glob function."""

    def test_none_returns_empty_list(self):
        """Should return empty list for None."""
        result = _read_file_list_or_glob(None, None)

        assert result == []

    def test_single_file_path(self, tmp_path: Path):
        """Should handle single file path."""
        file = tmp_path / "test.txt"
        file.touch()

        result = _read_file_list_or_glob(None, file)

        assert result == [file]

    def test_single_str_path(self, tmp_path: Path):
        """Should convert string to Path."""
        file = tmp_path / "test.txt"
        file.touch()

        result = _read_file_list_or_glob(None, str(file))

        assert result == [file]
        assert all(isinstance(p, Path) for p in result)

    def test_list_of_paths(self, tmp_path: Path):
        """Should handle list of paths."""
        files = [tmp_path / f"test{i}.txt" for i in range(3)]
        for f in files:
            f.touch()

        result = _read_file_list_or_glob(None, files)

        assert len(result) == 3
        assert all(p in result for p in files)

    def test_glob_pattern(self, tmp_path: Path):
        """Should expand glob patterns."""
        # Create multiple files
        for i in range(3):
            (tmp_path / f"data_{i}.txt").touch()

        pattern = str(tmp_path / "data_*.txt")

        result = _read_file_list_or_glob(None, pattern)

        assert len(result) == 3
        assert all(p.name.startswith("data_") for p in result)

    def test_glob_pattern_in_list(self, tmp_path: Path):
        """Should expand glob pattern provided as single-item list."""
        for i in range(3):
            (tmp_path / f"data_{i}.txt").touch()

        pattern = str(tmp_path / "data_*.txt")

        result = _read_file_list_or_glob(None, [pattern])

        assert len(result) == 3

    def test_text_file_with_list(self, tmp_path: Path):
        """Should read file list from text file."""
        # Create target files
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.touch()
        file2.touch()

        # Create list file
        list_file = tmp_path / "files.txt"
        list_file.write_text(f"{file1.name}\n{file2.name}\n")

        result = _read_file_list_or_glob(None, list_file)

        assert len(result) == 2
        assert file1 in result
        assert file2 in result

    def test_text_file_with_absolute_paths(self, tmp_path: Path):
        """Should handle absolute paths in text file."""
        file1 = tmp_path / "file1.txt"
        file2 = tmp_path / "file2.txt"
        file1.touch()
        file2.touch()

        list_file = tmp_path / "files.txt"
        list_file.write_text(f"{file1}\n{file2}\n")

        result = _read_file_list_or_glob(None, list_file)

        assert len(result) == 2
        assert file1 in result

    def test_text_file_ignores_empty_lines(self, tmp_path: Path):
        """Should skip empty lines in text file."""
        file1 = tmp_path / "file1.txt"
        file1.touch()

        list_file = tmp_path / "files.txt"
        list_file.write_text(f"{file1.name}\n\n  \n{file1.name}\n")

        result = _read_file_list_or_glob(None, list_file)

        # Should only have 2 entries even though file1 appears twice
        assert len(result) == 2

    def test_nonexistent_single_file(self, tmp_path: Path):
        """Should return path even if file doesn't exist."""
        nonexistent = tmp_path / "nonexistent.txt"

        result = _read_file_list_or_glob(None, nonexistent)

        # Should still return the path
        assert result == [nonexistent]

    def test_rejects_dict(self):
        """Should reject dict input."""
        with pytest.raises(TypeError, match="Expected string, Path, or list"):
            _read_file_list_or_glob(None, {"invalid": "dict"})


class TestConvertPathsToStrings:
    """Tests for convert_paths_to_strings function."""

    def test_single_path(self):
        """Should convert single Path to string."""
        path = Path("/some/path")

        result = convert_paths_to_strings(path)

        assert result == "/some/path"
        assert isinstance(result, str)

    def test_dict_with_paths(self):
        """Should convert Paths in dict."""
        data = {
            "file": Path("/file.txt"),
            "number": 42,
            "string": "text",
        }

        result = convert_paths_to_strings(data)

        assert result["file"] == "/file.txt"
        assert result["number"] == 42
        assert result["string"] == "text"

    def test_nested_dict(self):
        """Should recursively convert nested dicts."""
        data = {
            "outer": {
                "inner": Path("/nested/path"),
                "value": 123,
            }
        }

        result = convert_paths_to_strings(data)

        assert result["outer"]["inner"] == "/nested/path"
        assert result["outer"]["value"] == 123

    def test_list_with_paths(self):
        """Should convert Paths in list."""
        data = [
            Path("/path1"),
            "string",
            Path("/path2"),
            42,
        ]

        result = convert_paths_to_strings(data)

        assert result[0] == "/path1"
        assert result[1] == "string"
        assert result[2] == "/path2"
        assert result[3] == 42

    def test_dict_with_list_of_paths(self):
        """Should handle dict containing list of Paths."""
        data = {
            "files": [Path("/f1"), Path("/f2")],
            "count": 2,
        }

        result = convert_paths_to_strings(data)

        assert result["files"] == ["/f1", "/f2"]
        assert result["count"] == 2

    def test_preserves_non_path_types(self):
        """Should not modify non-Path types."""
        data = {
            "int": 42,
            "float": 3.14,
            "bool": True,
            "none": None,
            "str": "text",
        }

        result = convert_paths_to_strings(data)

        assert result == data


class TestFormatSummarySection:
    """Tests for format_summary_section function."""

    def test_basic_section(self):
        """Should format basic section."""
        items = {
            "name": "test",
            "value": 42,
        }

        lines = format_summary_section("Test Section", items)

        assert lines[0] == "Test Section"
        assert "=" in lines[1]  # Separator line
        assert "name: test" in " ".join(lines)
        assert "value: 42" in " ".join(lines)

    def test_custom_width(self):
        """Should respect custom width."""
        items = {"key": "value"}

        lines = format_summary_section("Title", items, max_width=40)

        # Separator should be 40 chars
        assert len(lines[1]) == 40

    def test_empty_items(self):
        """Should handle empty items dict."""
        items = {}

        lines = format_summary_section("Empty", items)

        # Should still have title and separator
        assert len(lines) >= 3
        assert lines[0] == "Empty"


class TestParseAlgorithmOverrides:
    """Tests for _parse_algorithm_overrides function."""

    def test_none_returns_empty(self):
        """Should return empty dict for None."""
        result = _parse_algorithm_overrides(None, 8882)

        assert result == {}

    def test_parses_frame_specific_overrides(self, tmp_path: Path):
        """Should extract frame-specific overrides."""
        overrides_file = tmp_path / "overrides.json"
        data = {
            "8882": {"downsample_factor": 5},
            "1234": {"downsample_factor": 10},
        }
        overrides_file.write_text(json.dumps(data))

        result = _parse_algorithm_overrides(overrides_file, 8882)

        assert result == {"downsample_factor": 5}

    def test_handles_data_wrapper(self, tmp_path: Path):
        """Should handle overrides wrapped in 'data' key."""
        overrides_file = tmp_path / "overrides.json"
        data = {
            "data": {
                "8882": {"downsample_factor": 5},
            }
        }
        overrides_file.write_text(json.dumps(data))

        result = _parse_algorithm_overrides(overrides_file, 8882)

        assert result == {"downsample_factor": 5}

    def test_returns_empty_for_missing_frame(self, tmp_path: Path):
        """Should return empty dict if frame not in overrides."""
        overrides_file = tmp_path / "overrides.json"
        data = {
            "1234": {"downsample_factor": 5},
        }
        overrides_file.write_text(json.dumps(data))

        result = _parse_algorithm_overrides(overrides_file, 9999)

        assert result == {}

    def test_handles_string_frame_id(self, tmp_path: Path):
        """Should convert frame_id to string for lookup."""
        overrides_file = tmp_path / "overrides.json"
        data = {
            "8882": {"downsample_factor": 5},
        }
        overrides_file.write_text(json.dumps(data))

        # Pass int frame_id
        result = _parse_algorithm_overrides(overrides_file, 8882)

        assert result == {"downsample_factor": 5}

    def test_complex_overrides(self, tmp_path: Path):
        """Should handle complex nested overrides."""
        overrides_file = tmp_path / "overrides.json"
        data = {
            "8882": {
                "downsample_factor": 5,
                "cal_method": "hanning_fft",
                "savitsky_goley": {
                    "window_x_size": 200,
                },
            }
        }
        overrides_file.write_text(json.dumps(data))

        result = _parse_algorithm_overrides(overrides_file, 8882)

        assert result["downsample_factor"] == 5
        assert result["cal_method"] == "hanning_fft"
        assert result["savitsky_goley"]["window_x_size"] == 200


class TestPathValidators:
    """Tests for path validator annotated types."""

    def test_required_path_accepts_path(self):
        """RequiredPath should accept Path."""
        from pydantic import BaseModel

        from cal_disp.config._utils import RequiredPath

        class TestModel(BaseModel):
            path: RequiredPath

        model = TestModel(path=Path("/test"))
        assert model.path == Path("/test")

    def test_required_path_accepts_string(self):
        """RequiredPath should convert string to Path."""
        from pydantic import BaseModel

        from cal_disp.config._utils import RequiredPath

        class TestModel(BaseModel):
            path: RequiredPath

        model = TestModel(path="/test")
        assert model.path == Path("/test")
        assert isinstance(model.path, Path)

    def test_required_path_rejects_none(self):
        """RequiredPath should reject None."""
        from pydantic import BaseModel, ValidationError

        from cal_disp.config._utils import RequiredPath

        class TestModel(BaseModel):
            path: RequiredPath

        with pytest.raises(ValidationError):
            TestModel(path=None)

    def test_optional_path_accepts_none(self):
        """OptionalPath should accept None."""
        from pydantic import BaseModel

        from cal_disp.config._utils import OptionalPath

        class TestModel(BaseModel):
            path: OptionalPath = None

        model = TestModel()
        assert model.path is None

    def test_optional_path_accepts_path(self):
        """OptionalPath should accept Path."""
        from pydantic import BaseModel

        from cal_disp.config._utils import OptionalPath

        class TestModel(BaseModel):
            path: OptionalPath = None

        model = TestModel(path=Path("/test"))
        assert model.path == Path("/test")

    def test_directory_path_accepts_empty_string(self):
        """DirectoryPath should accept empty string as current dir."""
        from pydantic import BaseModel

        from cal_disp.config._utils import DirectoryPath

        class TestModel(BaseModel):
            path: DirectoryPath

        model = TestModel(path="")
        assert model.path == Path()
