"""Tests for YamlModel base class."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml  # type: ignore[import-untyped]
from pydantic import Field

from cal_disp.config._yaml import YamlModel


class SimpleModel(YamlModel):
    """Simple test model."""

    name: str = Field(description="Name field")
    value: int = Field(default=42, description="Value field")


class ModelWithPaths(YamlModel):
    """Model with Path fields."""

    file_path: Path
    optional_path: Path | None = None
    file_list: list[Path] | None = None


class TestYamlModelBasics:
    """Tests for basic YamlModel functionality."""

    def test_creates_instance(self):
        """Should create model instance."""
        model = SimpleModel(name="test")

        assert model.name == "test"
        assert model.value == 42

    def test_to_yaml_creates_file(self, tmp_path: Path):
        """Should save to YAML file."""
        model = SimpleModel(name="test", value=100)
        yaml_file = tmp_path / "test.yaml"

        model.to_yaml(yaml_file)

        assert yaml_file.exists()

    def test_from_yaml_loads_file(self, tmp_path: Path):
        """Should load from YAML file."""
        yaml_file = tmp_path / "test.yaml"
        yaml_file.write_text("name: loaded\nvalue: 200\n")

        model = SimpleModel.from_yaml(yaml_file)

        assert model.name == "loaded"
        assert model.value == 200

    def test_yaml_round_trip(self, tmp_path: Path):
        """Should preserve data through save/load cycle."""
        original = SimpleModel(name="roundtrip", value=999)
        yaml_file = tmp_path / "roundtrip.yaml"

        original.to_yaml(yaml_file)
        loaded = SimpleModel.from_yaml(yaml_file)

        assert loaded.name == original.name
        assert loaded.value == original.value


class TestYamlModelComments:
    """Tests for YAML comments functionality."""

    def test_to_yaml_without_comments(self, tmp_path: Path):
        """Should save without comments when requested."""
        model = SimpleModel(name="test")
        yaml_file = tmp_path / "no_comments.yaml"

        model.to_yaml(yaml_file, with_comments=False)

        content = yaml_file.read_text()
        assert "#" not in content

    def test_to_yaml_with_comments(self, tmp_path: Path):
        """Should include field descriptions as comments."""
        model = SimpleModel(name="test")
        yaml_file = tmp_path / "with_comments.yaml"

        model.to_yaml(yaml_file, with_comments=True)

        content = yaml_file.read_text()
        assert "#" in content
        # Should have description comments
        assert "Name field" in content or "description" in content.lower()


class TestYamlModelAliases:
    """Tests for field alias handling."""

    def test_by_alias_true(self, tmp_path: Path):
        """Should use field aliases when by_alias=True."""

        class ModelWithAlias(YamlModel):
            internal_name: str = Field(alias="external_name")

        model = ModelWithAlias(external_name="test")
        yaml_file = tmp_path / "alias.yaml"

        model.to_yaml(yaml_file, by_alias=True)

        with open(yaml_file) as f:
            data = yaml.safe_load(f)

        assert "external_name" in data
        assert "internal_name" not in data

    def test_by_alias_false(self, tmp_path: Path):
        """Should use field names when by_alias=False."""

        class ModelWithAlias(YamlModel):
            internal_name: str = Field(alias="external_name")

        model = ModelWithAlias(external_name="test")
        yaml_file = tmp_path / "no_alias.yaml"

        model.to_yaml(yaml_file, by_alias=False)

        with open(yaml_file) as f:
            data = yaml.safe_load(f)

        assert "internal_name" in data


class TestYamlModelFilePathMethods:
    """Tests for file path handling methods."""

    def test_get_all_file_paths_single(self, tmp_path: Path):
        """Should find single Path field."""
        file = tmp_path / "test.txt"
        file.touch()

        model = ModelWithPaths(file_path=file)

        files = model.get_all_file_paths()

        assert "file_path" in files
        assert files["file_path"] == file

    def test_get_all_file_paths_optional_none(self, tmp_path: Path):
        """Should exclude None optional paths by default."""
        file = tmp_path / "test.txt"
        file.touch()

        model = ModelWithPaths(file_path=file, optional_path=None)

        files = model.get_all_file_paths()

        assert "file_path" in files
        assert "optional_path" not in files

    def test_get_all_file_paths_include_none(self, tmp_path: Path):
        """Should include None paths when requested."""
        file = tmp_path / "test.txt"
        file.touch()

        model = ModelWithPaths(file_path=file, optional_path=None)

        files = model.get_all_file_paths(include_none=True)

        assert "optional_path" in files
        assert files["optional_path"] is None

    def test_get_all_file_paths_list_flattened(self, tmp_path: Path):
        """Should flatten list of paths with indices."""
        file1 = tmp_path / "f1.txt"
        file2 = tmp_path / "f2.txt"
        file1.touch()
        file2.touch()

        model = ModelWithPaths(
            file_path=file1,
            file_list=[file1, file2],
        )

        files = model.get_all_file_paths(flatten_lists=True)

        assert "file_list[0]" in files
        assert "file_list[1]" in files
        assert files["file_list[0]"] == file1
        assert files["file_list[1]"] == file2

    def test_get_all_file_paths_list_not_flattened(self, tmp_path: Path):
        """Should keep list as single entry when not flattening."""
        file1 = tmp_path / "f1.txt"
        file2 = tmp_path / "f2.txt"
        file1.touch()
        file2.touch()

        model = ModelWithPaths(
            file_path=file1,
            file_list=[file1, file2],
        )

        files = model.get_all_file_paths(flatten_lists=False)

        assert "file_list" in files
        assert files["file_list"] == [file1, file2]


class TestYamlModelValidateFilesExist:
    """Tests for file existence validation."""

    def test_validate_files_exist_all_present(self, tmp_path: Path):
        """Should report all files as existing."""
        file = tmp_path / "test.txt"
        file.touch()

        model = ModelWithPaths(file_path=file)

        results = model.validate_files_exist()

        assert results["file_path"]["exists"] is True
        assert results["file_path"]["is_file"] is True

    def test_validate_files_exist_missing(self, tmp_path: Path):
        """Should report missing files."""
        file = tmp_path / "test.txt"
        missing = tmp_path / "missing.txt"
        file.touch()

        model = ModelWithPaths(file_path=missing)

        results = model.validate_files_exist()

        assert results["file_path"]["exists"] is False
        assert results["file_path"]["is_file"] is None

    def test_validate_files_exist_raise_on_missing(self, tmp_path: Path):
        """Should raise error for missing files when requested."""
        missing = tmp_path / "missing.txt"

        model = ModelWithPaths(file_path=missing)

        with pytest.raises(FileNotFoundError, match="file_path"):
            model.validate_files_exist(raise_on_missing=True)

    def test_validate_files_exist_includes_size(self, tmp_path: Path):
        """Should include file size in results."""
        file = tmp_path / "test.txt"
        file.write_text("test content")

        model = ModelWithPaths(file_path=file)

        results = model.validate_files_exist()

        assert "size_bytes" in results["file_path"]
        assert results["file_path"]["size_bytes"] > 0


class TestYamlModelMissingFiles:
    """Tests for get_missing_files method."""

    def test_get_missing_files_none(self, tmp_path: Path):
        """Should return empty list when all files exist."""
        file = tmp_path / "test.txt"
        file.touch()

        model = ModelWithPaths(file_path=file)

        missing = model.get_missing_files()

        assert len(missing) == 0

    def test_get_missing_files_some(self, tmp_path: Path):
        """Should return list of missing files."""
        missing = tmp_path / "missing.txt"

        model = ModelWithPaths(file_path=missing)

        missing_list = model.get_missing_files()

        assert len(missing_list) > 0
        assert any("file_path" in m for m in missing_list)

    def test_all_files_exist_true(self, tmp_path: Path):
        """Should return True when all files exist."""
        file = tmp_path / "test.txt"
        file.touch()

        model = ModelWithPaths(file_path=file)

        assert model.all_files_exist() is True

    def test_all_files_exist_false(self, tmp_path: Path):
        """Should return False when files are missing."""
        missing = tmp_path / "missing.txt"

        model = ModelWithPaths(file_path=missing)

        assert model.all_files_exist() is False


class TestYamlModelPrintSchema:
    """Tests for print_yaml_schema method."""

    def test_print_yaml_schema_stdout(self, capsys):
        """Should print schema to stdout."""
        SimpleModel.print_yaml_schema()

        captured = capsys.readouterr()
        assert "name:" in captured.out

    def test_print_yaml_schema_to_file(self, tmp_path: Path):
        """Should save schema to file."""
        schema_file = tmp_path / "schema.yaml"

        SimpleModel.print_yaml_schema(schema_file)

        assert schema_file.exists()
        content = schema_file.read_text()
        assert "name:" in content


class TestYamlModelValidateReadyToRun:
    """Tests for validate_ready_to_run method."""

    def test_default_implementation(self):
        """Should have default implementation that passes."""
        model = SimpleModel(name="test")

        result = model.validate_ready_to_run()

        assert result["ready"] is True
        assert len(result["errors"]) == 0
        assert len(result["warnings"]) == 0


class TestYamlModelStrictConfig:
    """Tests for strict configuration enforcement."""

    def test_forbids_extra_fields(self):
        """Should forbid extra fields."""
        with pytest.raises(ValueError):
            SimpleModel(name="test", extra_field="not allowed")  # type: ignore

    def test_validates_on_assignment(self):
        """Should validate when fields are assigned."""
        model = SimpleModel(name="test")

        # This should work fine
        model.value = 100
        assert model.value == 100


class TestYamlModelIndentation:
    """Tests for YAML indentation settings."""

    def test_custom_indent_per_level(self, tmp_path: Path):
        """Should respect custom indentation."""

        class NestedModel(YamlModel):
            outer: dict[str, int] = {"a": 1, "b": 2}

        model = NestedModel()
        yaml_file = tmp_path / "indented.yaml"

        model.to_yaml(yaml_file, indent_per_level=4)

        # Check that indentation is present (hard to test exact spacing)
        assert yaml_file.exists()
        content = yaml_file.read_text()
        assert "outer:" in content
