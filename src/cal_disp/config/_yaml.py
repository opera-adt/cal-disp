from __future__ import annotations

import json
import sys
import textwrap
from io import StringIO
from os import PathLike
from pathlib import Path
from typing import Any, Dict, List, Optional, TextIO, TypedDict, Union, cast

from pydantic import BaseModel, ConfigDict
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap

Filename = Union[str, PathLike]

# Standard configs
STRICT_CONFIG = ConfigDict(
    extra="forbid",
    validate_assignment=True,
)

STRICT_CONFIG_WITH_ALIASES = ConfigDict(
    extra="forbid",
    validate_assignment=True,
    populate_by_name=True,
)


class ValidationResult(TypedDict):
    """Standard validation result structure."""

    # Using dict subclass instead of TypedDict for runtime flexibility
    ready: bool
    errors: List[str]
    warnings: List[str]


class YamlModel(BaseModel):
    """Pydantic model that can be exported to yaml."""

    model_config = STRICT_CONFIG

    def to_yaml(
        self,
        output_path: Union[Filename, TextIO],
        with_comments: bool = True,
        by_alias: bool = True,
        indent_per_level: int = 2,
    ):
        """Save configuration as a yaml file.

        Used to record the default-filled version of a supplied yaml.

        Parameters
        ----------
        output_path : Pathlike
            Path to the yaml file to save.
        with_comments : bool, default = False.
            Whether to add comments containing the type/descriptions to all fields.
        by_alias : bool, default = False.
            Whether to use the alias names for the fields.
            Passed to pydantic's ``to_json`` method.
            https://docs.pydantic.dev/usage/exporting_models/#modeljson
        indent_per_level : int, default = 2
            Number of spaces to indent per level.

        """
        yaml_obj = self._to_yaml_obj(by_alias=by_alias)

        if with_comments:
            _add_comments(
                yaml_obj,
                self.model_json_schema(by_alias=by_alias),
                indent_per_level=indent_per_level,
            )

        y = YAML()
        # https://yaml.readthedocs.io/en/latest/detail.html#indentation-of-block-sequences
        y.indent(
            offset=indent_per_level,
            mapping=indent_per_level,
            # It is best to always have sequence >= offset + 2 but this is not enforced
            # not following this advice might lead to invalid output.
            sequence=indent_per_level + 2,
        )
        if hasattr(output_path, "write"):
            y.dump(yaml_obj, output_path)
        else:
            with open(output_path, "w") as f:
                y.dump(yaml_obj, f)

    @classmethod
    def from_yaml(cls, yaml_path: Filename):
        """Load a configuration from a yaml file.

        Parameters
        ----------
        yaml_path : Pathlike
            Path to the yaml file to load.

        Returns
        -------
        Config
            Workflow configuration

        """
        y = YAML(typ="safe")
        with open(yaml_path) as f:
            data = y.load(f)

        return cls(**data)

    @classmethod
    def print_yaml_schema(
        cls,
        output_path: Union[Filename, TextIO] = sys.stdout,
        indent_per_level: int = 2,
    ):
        """Print/save an empty configuration with defaults filled in.

        Ignores the required `input_file_list` input, so a user can
        inspect all fields.

        Parameters
        ----------
        output_path : Pathlike
            Path or stream to save to the yaml file to.
            By default, prints to stdout.
        indent_per_level : int, default = 2
            Number of spaces to indent per level.

        """
        cls.model_construct().to_yaml(
            output_path, with_comments=True, indent_per_level=indent_per_level
        )

    def _to_yaml_obj(self, by_alias: bool = True) -> CommentedMap:
        # Make the YAML object to add comments to
        # We can't just do `dumps` for some reason, need a stream
        y = YAML()
        ss = StringIO()
        y.dump(json.loads(self.model_dump_json(by_alias=by_alias)), ss)
        return y.load(ss.getvalue())

    def get_all_file_paths(
        self, include_none: bool = False, flatten_lists: bool = True
    ) -> Dict[str, Union[Path, List[Path]]]:
        """Get all Path fields from the model.

        Parameters
        ----------
        include_none : bool, default=False
            Include fields with None values.
        flatten_lists : bool, default=True
            Flatten list fields to individual entries with indices.

        Returns
        -------
        dict
            Mapping of field names to Path objects.

        """
        files: Dict[str, Path | list[Path]] = {}

        for field_name, field_info in self.model_fields.items():
            value = getattr(self, field_name)

            # Skip None values if requested
            if value is None and not include_none:
                continue

            # Check if field is Path or Optional[Path]
            if self._is_path_field(field_info):
                if value is not None:
                    files[field_name] = value

            # Check if field is List[Path]
            elif self._is_path_list_field(field_info):
                if value:
                    if flatten_lists:
                        for i, path in enumerate(value):
                            files[f"{field_name}[{i}]"] = path
                    else:
                        files[field_name] = value

        return files

    @staticmethod
    def _is_path_field(field_info) -> bool:
        """Check if field is a Path type."""
        from pathlib import Path
        from typing import Union, get_args, get_origin

        annotation = field_info.annotation

        # Handle Annotated[Path, ...]
        if get_origin(annotation) is not None:
            # Check if it's Annotated
            if hasattr(annotation, "__metadata__"):
                # Get the actual type from Annotated
                args = get_args(annotation)
                if args:
                    annotation = args[0]

        # Direct Path type
        if annotation is Path:
            return True

        # Optional[Path] or Union[Path, None]
        origin = get_origin(annotation)
        if origin is Union:
            args = get_args(annotation)
            return Path in args or any(arg is Path for arg in args)

        return False

    @staticmethod
    def _is_path_list_field(field_info) -> bool:
        """Check if field is a List[Path] type."""
        from pathlib import Path
        from typing import Union, get_args, get_origin

        annotation = field_info.annotation
        origin = get_origin(annotation)

        # Check if it's Optional[List[Path]]
        if origin is Union:
            args = get_args(annotation)
            for arg in args:
                if arg is type(None):
                    continue
                if get_origin(arg) in (list, List):
                    list_args = get_args(arg)
                    if list_args:
                        first_arg = list_args[0]
                        if first_arg is Path:
                            return True

        # Check if it's List[Path]
        if origin in (list, List):
            args = get_args(annotation)
            if not args:
                return False
            first_arg = args[0]
            is_path_type: bool = first_arg is Path
            return is_path_type

        return False

    def validate_files_exist(
        self, raise_on_missing: bool = False
    ) -> Dict[str, Dict[str, Any]]:
        """Validate all file paths exist on disk.

        Parameters
        ----------
        raise_on_missing : bool, default=False
            If True, raise FileNotFoundError on first missing file.

        Returns
        -------
        dict
            Detailed status for each file including existence, size, etc.

        Raises
        ------
        FileNotFoundError
            If raise_on_missing=True and any file is missing.

        """
        results: Dict[str, Dict[str, Any]] = {}

        # flatten_lists=True guarantees Dict[str, Path]
        file_paths: Dict[str, Path] = cast(
            Dict[str, Path], self.get_all_file_paths(flatten_lists=True)
        )

        for field_name, path in file_paths.items():
            exists = path.exists()

            if not exists and raise_on_missing:
                raise FileNotFoundError(
                    f"Required file not found: {field_name} = {path}"
                )

            results[field_name] = {
                "exists": exists,
                "is_file": path.is_file() if exists else None,
                "is_dir": path.is_dir() if exists else None,
                "size_bytes": path.stat().st_size if exists else None,
                "absolute_path": str(path.absolute()),
            }

        return results

    def validate_ready_to_run(self) -> ValidationResult:
        """Check if configuration is ready to run."""
        return ValidationResult(ready=True, errors=[], warnings=[])

    def get_missing_files(self) -> List[str]:
        """Get list of missing file paths."""
        return [
            f"{name}: {info['absolute_path']}"
            for name, info in self.validate_files_exist().items()
            if not info["exists"]
        ]

    def all_files_exist(self) -> bool:
        """Check if all files exist."""
        return len(self.get_missing_files()) == 0


def _add_comments(
    loaded_yaml: CommentedMap,
    schema: Dict[str, Any],
    indent: int = 0,
    definitions: Optional[dict] = None,
    # variable specifying how much to indent per level
    indent_per_level: int = 2,
):
    """Add comments above each YAML field using the pydantic model schema."""
    # Definitions are in schemas that contain nested pydantic Models
    defs = schema.get("$defs") if definitions is None else definitions

    for key, val in schema["properties"].items():
        reference = ""
        # Get sub-schema if it exists
        if "$ref" in val:
            # At top level, example is 'outputs': {'$ref': '#/defs/Outputs'}
            reference = val["$ref"]
        elif "allOf" in val:
            # within 'defs', it looks like
            #  'allOf': [{'$ref': '#/defs/HalfWindow'}]
            reference = val["allOf"][0]["$ref"]

        ref_key = reference.split("/")[-1]
        if ref_key:  # The current property is a reference to something else
            if "enum" in defs[ref_key]:  # type: ignore[index]
                # This is just an Enum, not a sub schema.
                # Overwrite the value with the referenced value
                val = defs[ref_key]  # type: ignore[index] # noqa: PLW2901
            else:
                # The reference is a sub schema, so we need to recurse
                sub_schema = defs[ref_key]  # type: ignore[index]
                # Get the sub-model
                sub_loaded_yaml = loaded_yaml[key]
                # recurse on the sub-model
                _add_comments(
                    sub_loaded_yaml,
                    sub_schema,
                    indent=indent + indent_per_level,
                    definitions=defs,
                    indent_per_level=indent_per_level,
                )
                continue

        # add each description along with the type information
        desc = val.get("description", "No description.")
        desc = "\n".join(
            textwrap.wrap(
                f"{desc}.",
                width=90,
                subsequent_indent=" " * indent_per_level,
            )
        )
        if "anyOf" in val:
            #   'anyOf': [{'type': 'string'}, {'type': 'null'}],
            # Join the options with a pipe, like Python types
            type_str = " | ".join(t.get("type", "None") for t in val["anyOf"]).replace(
                "null", "None"
            )
        elif "const" in val:
            type_str = val["const"]
        else:
            type_str = val.get("type", "Any")
        type_line = f"\n  Type: {type_str}."
        choices = f"\n  Options: {val['enum']}." if "enum" in val else ""

        # Combine the description/type/choices as the YAML comment
        comment = f"{desc}{type_line}{choices}".replace("..", ".")

        # Prepend the required label for fields that are required
        is_required = key in schema.get("required", [])
        if is_required:
            comment = "REQUIRED: " + comment

        # This method comes from here
        # https://yaml.readthedocs.io/en/latest/detail.html#round-trip-including-comments
        loaded_yaml.yaml_set_comment_before_after_key(
            key,
            comment,
            indent=indent,
        )
