from __future__ import annotations

import json
import sys
import textwrap
from io import StringIO
from os import PathLike
from pathlib import Path
from typing import (
    Any,
    ClassVar,
    Dict,
    List,
    Optional,
    TextIO,
    TypedDict,
    Union,
    cast,
    get_args,
    get_origin,
)

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

    ready: bool
    errors: List[str]
    warnings: List[str]


def _get_placeholder_value(annotation):
    """Get a placeholder value for a type."""
    origin = get_origin(annotation)

    # Handle Annotated types
    if hasattr(annotation, "__metadata__"):
        args = get_args(annotation)
        if args:
            annotation = args[0]
            origin = get_origin(annotation)

    if annotation is str:
        return ""
    elif annotation is int:
        return 0
    elif annotation is bool:
        return False
    elif annotation is Path:
        return Path()
    elif origin in (list, List):
        return []
    elif origin in (dict, Dict):
        return {}
    else:
        return None


class YamlModel(BaseModel):
    """Pydantic model that can be exported to yaml."""

    _yaml_root_key: ClassVar[Optional[str]] = None
    model_config = STRICT_CONFIG

    def to_yaml(
        self,
        output_path: Union[Filename, TextIO],
        with_comments: bool = True,
        by_alias: bool = True,
        indent_per_level: int = 2,
    ):
        """Save configuration as a yaml file with support for wrapped comments."""
        yaml_obj = self._to_yaml_obj(by_alias=by_alias)

        # Get the root key safely from the ClassVar
        root_key = getattr(self, "_yaml_root_key", None)

        if with_comments:
            schema = self.model_json_schema(by_alias=by_alias)

            # COMMENT LOGIC:
            # If wrapped, we point the comment logic at the values inside the key
            # and increase the base indent.
            target_yaml = yaml_obj
            comment_indent = 0

            if root_key and root_key in yaml_obj:
                target_yaml = yaml_obj[root_key]
                comment_indent = indent_per_level

                # Add a high-level comment to the top-level key itself
                desc = (
                    self.__doc__.split("\n")[0].strip()
                    if self.__doc__
                    else "Workflow configuration."
                )
                yaml_obj.yaml_set_comment_before_after_key(
                    root_key, f"REQUIRED: {desc}", indent=0
                )

            # Recursive call to add comments for all fields
            _add_comments(
                target_yaml,
                schema,
                indent=comment_indent,
                definitions=schema.get("$defs"),
                indent_per_level=indent_per_level,
            )

        y = YAML()
        y.indent(
            offset=indent_per_level,
            mapping=indent_per_level,
            sequence=indent_per_level + 2,
        )
        if hasattr(output_path, "write"):
            y.dump(yaml_obj, output_path)
        else:
            with open(output_path, "w") as f:
                y.dump(yaml_obj, f)

    @classmethod
    def from_yaml(cls, yaml_path: Filename):
        """Load YAML, unwrapping the top-level key if it matches cls.name."""
        y = YAML(typ="safe")
        with open(yaml_path) as f:
            data = y.load(f)

        # Check the unique attribute name
        root_key = getattr(cls, "_yaml_root_key", None)
        if root_key and isinstance(data, dict) and root_key in data:
            data = data[root_key]

        return cls(**data)

    @classmethod
    def print_yaml_schema(
        cls,
        output_path: Union[Filename, TextIO, None] = None,
        indent_per_level: int = 2,
    ):
        """Print/save an empty configuration with defaults filled in.

        Shows all fields including required ones (with empty placeholders).
        """
        # Resolve stdout at runtime so pytest can capture it
        if output_path is None:
            output_path = sys.stdout

        # Get placeholder values for required fields
        field_values = {}
        for field_name, field_info in cls.model_fields.items():
            if field_info.is_required():
                field_values[field_name] = _get_placeholder_value(field_info.annotation)

        cls.model_construct(**field_values).to_yaml(
            output_path, with_comments=True, indent_per_level=indent_per_level
        )

    def _to_yaml_obj(self, by_alias: bool = True) -> CommentedMap:
        """Convert model to YAML, wrapping in 'name' if defined."""
        y = YAML()
        ss = StringIO()

        # Get the flat dictionary
        json_data = self.model_dump_json(
            by_alias=by_alias,
            exclude_unset=False,
            exclude_none=False,
        )
        data = json.loads(json_data)

        root_key = getattr(self, "_yaml_root_key", None)
        if root_key:
            data = {root_key: data}

        y.dump(data, ss)
        return y.load(ss.getvalue())

    def get_all_file_paths(
        self, include_none: bool = False, flatten_lists: bool = True
    ) -> Dict[str, Union[Path, List[Path]]]:
        """Get all Path fields from the model."""
        files: Dict[str, Path | list[Path]] = {}

        for field_name, field_info in self.__class__.model_fields.items():
            value = getattr(self, field_name)

            # Skip None values if requested
            if value is None and not include_none:
                continue

            # Check if field is Path or Optional[Path]
            if self._is_path_field(field_info):
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
        import types
        from typing import Union

        annotation = field_info.annotation

        # Handle Annotated[Path, ...]
        origin = get_origin(annotation)
        if origin is not None and hasattr(annotation, "__metadata__"):
            args = get_args(annotation)
            if args:
                annotation = args[0]

        # Direct Path type
        if annotation is Path:
            return True

        # Optional[Path] or Path | None
        origin = get_origin(annotation)
        if origin is Union or (
            hasattr(types, "UnionType") and origin is types.UnionType
        ):
            args = get_args(annotation)
            return Path in args

        return False

    @staticmethod
    def _is_path_list_field(field_info) -> bool:
        """Check if field is a List[Path] type."""
        import types
        from typing import Union

        annotation = field_info.annotation
        origin = get_origin(annotation)

        # Check if it's Optional[List[Path]] or List[Path] | None
        if origin is Union or (
            hasattr(types, "UnionType") and origin is types.UnionType
        ):
            args = get_args(annotation)
            for arg in args:
                if arg is type(None):
                    continue
                if get_origin(arg) in (list, List):
                    list_args = get_args(arg)
                    if list_args and list_args[0] is Path:
                        return True

        # Check if it's List[Path]
        if origin in (list, List):
            args = get_args(annotation)
            if args and args[0] is Path:
                return True

        return False

    def validate_files_exist(
        self, raise_on_missing: bool = False
    ) -> Dict[str, Dict[str, Any]]:
        """Validate all file paths exist on disk."""
        results: Dict[str, Dict[str, Any]] = {}

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
    loaded_yaml: Any,
    schema: Dict[str, Any],
    indent: int = 0,
    definitions: Optional[dict] = None,
    indent_per_level: int = 2,
):
    """Add comments above each YAML field using the pydantic model schema."""
    if definitions is None:
        definitions = schema.get("$defs", {})

    properties = schema.get("properties", {})

    for key, val in properties.items():
        # 1. Resolve references and find the actual "value" schema
        current_val = val
        reference = None

        # Unwrap anyOf/allOf (common in Optional fields)
        if "anyOf" in current_val:
            # Find the branch that isn't 'null' and has a reference
            for branch in current_val["anyOf"]:
                if "$ref" in branch:
                    reference = branch["$ref"]
                    break
        elif "$ref" in current_val:
            reference = current_val["$ref"]
        elif "allOf" in current_val:
            reference = current_val["allOf"][0].get("$ref")

        # 2. Handle Recursion if it's a nested model
        if reference:
            ref_key = reference.split("/")[-1]
            if ref_key in definitions:
                sub_schema = definitions[ref_key]
                # If the YAML has this key and it's a dict, dive in
                if key in loaded_yaml and isinstance(loaded_yaml[key], dict):
                    _add_comments(
                        loaded_yaml[key],
                        sub_schema,
                        indent=indent + indent_per_level,
                        definitions=definitions,
                        indent_per_level=indent_per_level,
                    )

        # 3. Process description for the current key
        desc = val.get("description", "No description.")
        # Wrap text to maintain clean formatting in YAML
        wrapped_desc = "\n".join(
            textwrap.wrap(
                f"{desc}.",
                width=90,
                subsequent_indent=" " * (indent + indent_per_level),
            )
        )

        # Determine type string for the comment
        if "anyOf" in val:
            type_str = " | ".join(
                t.get("type", t.get("$ref", "").split("/")[-1]) for t in val["anyOf"]
            ).replace("NoneType", "null")
        else:
            type_str = val.get("type", "Object")

        comment = f"{wrapped_desc}\n  Type: {type_str}."
        if key in schema.get("required", []):
            comment = "REQUIRED: " + comment

        # 4. Set the comment
        if hasattr(loaded_yaml, "yaml_set_comment_before_after_key"):
            try:
                loaded_yaml.yaml_set_comment_before_after_key(
                    key, comment, indent=indent
                )
            except (KeyError, TypeError):
                pass
