"""Tests for PGE run configuration models."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml  # type: ignore[import-untyped]

from cal_disp.config.pge_runconfig import (
    OutputOptions,
    PrimaryExecutable,
    ProductPathGroup,
    RunConfig,
)


class TestPrimaryExecutable:
    """Tests for PrimaryExecutable configuration."""

    def test_default_product_type(self):
        """Should have CAL_DISP as default product type."""
        exec_config = PrimaryExecutable()

        assert exec_config.product_type == "CAL_DISP"

    def test_custom_product_type(self):
        """Should accept custom product type."""
        exec_config = PrimaryExecutable(product_type="CUSTOM")

        assert exec_config.product_type == "CUSTOM"

    def test_forbids_extra_fields(self):
        """Should forbid extra fields."""
        with pytest.raises(ValueError):
            PrimaryExecutable(extra_field="not allowed")  # type: ignore


class TestOutputOptions:
    """Tests for OutputOptions configuration."""

    def test_defaults(self):
        """Should have expected defaults."""
        options = OutputOptions()

        assert options.product_version == "1.0"
        assert options.output_format == "netcdf"
        assert options.compression is True

    def test_custom_version(self):
        """Should accept custom version."""
        options = OutputOptions(product_version="2.5")

        assert options.product_version == "2.5"

    def test_custom_format(self):
        """Should accept custom output format."""
        options = OutputOptions(output_format="hdf5")

        assert options.output_format == "hdf5"

    def test_disable_compression(self):
        """Should allow disabling compression."""
        options = OutputOptions(compression=False)

        assert options.compression is False


class TestProductPathGroup:
    """Tests for ProductPathGroup configuration."""

    def test_defaults(self):
        """Should have expected default paths."""
        paths = ProductPathGroup()

        assert paths.product_path == Path()
        assert paths.scratch_path == Path("./scratch")
        assert paths.output_path == Path("./output")

    def test_custom_paths(self, tmp_path: Path):
        """Should accept custom paths."""
        product = tmp_path / "product"
        scratch = tmp_path / "scratch"
        output = tmp_path / "output"

        paths = ProductPathGroup(
            product_path=product,
            scratch_path=scratch,
            output_path=output,
        )

        assert paths.product_path == product
        assert paths.scratch_path == scratch
        assert paths.output_path == output

    def test_sas_output_path_alias(self, tmp_path: Path):
        """Should accept sas_output_path as alias."""
        output = tmp_path / "output"

        paths = ProductPathGroup(sas_output_path=output)

        assert paths.output_path == output


class TestRunConfig:
    """Tests for complete RunConfig."""

    def test_requires_input_file_group(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should require InputFileGroup."""
        from cal_disp.config import InputFileGroup

        with pytest.raises(ValueError):
            RunConfig()  # Missing required input_file_group  # type: ignore

        # Should work with input_file_group
        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )
        )
        assert config.input_file_group is not None

    def test_creates_with_defaults(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should create with default values for optional fields."""
        from cal_disp.config import InputFileGroup

        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )
        )

        assert config.output_options.product_version == "1.0"
        assert config.primary_executable.product_type == "CAL_DISP"
        assert config.worker_settings.n_workers == 4

    def test_with_all_groups(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        tmp_path: Path,
    ):
        """Should accept all optional file groups."""
        from cal_disp.config import (
            DynamicAncillaryFileGroup,
            InputFileGroup,
            StaticAncillaryFileGroup,
        )

        los_file, dem_file = sample_static_layers

        override_file = tmp_path / "overrides.json"
        override_file.write_text("{}")

        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            dynamic_ancillary_group=DynamicAncillaryFileGroup(
                algorithm_parameters_file=sample_algorithm_params,
                static_los_file=los_file,
                static_dem_file=dem_file,
            ),
            static_ancillary_group=StaticAncillaryFileGroup(
                algorithm_parameters_overrides_json=override_file
            ),
        )

        assert config.dynamic_ancillary_group is not None
        assert config.static_ancillary_group is not None

    def test_create_directories(
        self,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should create all necessary directories."""
        from cal_disp.config import InputFileGroup

        product_dir = tmp_path / "product"
        scratch_dir = tmp_path / "scratch"
        output_dir = tmp_path / "output"

        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            product_path_group=ProductPathGroup(
                product_path=product_dir,
                scratch_path=scratch_dir,
                output_path=output_dir,
            ),
        )

        assert not product_dir.exists()
        assert not scratch_dir.exists()
        assert not output_dir.exists()

        config.create_directories()

        assert product_dir.exists()
        assert scratch_dir.exists()
        assert output_dir.exists()
        assert (scratch_dir / "tmp").exists()

    def test_validate_ready_to_run_missing_inputs(self, tmp_path: Path):
        """Should report errors for missing required inputs."""
        from cal_disp.config import InputFileGroup

        fake_grid_dir = tmp_path / "grid"
        fake_grid_dir.mkdir()  # create the directory

        # create a dummy .tenv8 file
        (fake_grid_dir / "dummy.tenv8").write_text("fake content")

        # Create with minimal inputs
        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=Path("fake.nc"),
                unr_grid_latlon_file=Path("grid_latlon_lookup_v0.2.txt"),
                unr_timeseries_dir=fake_grid_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )
        )

        validation = config.validate_ready_to_run()

        # Should have warnings about missing dynamic ancillary
        assert not validation["ready"] or len(validation["warnings"]) > 0

    def test_validate_ready_to_run_complete(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should pass validation with all required inputs."""
        from cal_disp.config import DynamicAncillaryFileGroup, InputFileGroup

        los_file, dem_file = sample_static_layers

        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            dynamic_ancillary_group=DynamicAncillaryFileGroup(
                algorithm_parameters_file=sample_algorithm_params,
                static_los_file=los_file,
                static_dem_file=dem_file,
            ),
        )

        validation = config.validate_ready_to_run()

        assert validation["ready"]
        assert len(validation["errors"]) == 0

    def test_summary(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should generate readable summary."""
        from cal_disp.config import InputFileGroup

        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )
        )

        summary = config.summary()

        assert "PGE Run Configuration" in summary
        assert "8882" in summary  # Frame ID
        assert "CAL_DISP" in summary  # Product type

    def test_to_workflow(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should convert to CalibrationWorkflow."""
        from cal_disp.config import DynamicAncillaryFileGroup, InputFileGroup

        los_file, dem_file = sample_static_layers

        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            dynamic_ancillary_group=DynamicAncillaryFileGroup(
                algorithm_parameters_file=sample_algorithm_params,
                static_los_file=los_file,
                static_dem_file=dem_file,
            ),
        )

        workflow = config.to_workflow()

        from cal_disp.config.workflow import CalibrationWorkflow

        assert isinstance(workflow, CalibrationWorkflow)
        assert workflow.input_options is not None
        assert workflow.input_options.frame_id == 8882

    def test_yaml_round_trip(
        self,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should save and load via YAML with wrapper key."""
        from cal_disp.config import InputFileGroup

        yaml_file = tmp_path / "runconfig.yaml"

        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )
        )

        # Save
        config.to_yaml(yaml_file)

        # Check wrapper key exists
        with open(yaml_file) as f:
            data = yaml.safe_load(f)
        assert "cal_disp_workflow" in data

        # Load
        loaded = RunConfig.from_yaml_file(yaml_file)

        assert loaded.input_file_group.frame_id == 8882

    def test_create_example(self):
        """Should create example configuration."""
        config = RunConfig.create_example()

        assert config.input_file_group.frame_id == 8882
        assert config.dynamic_ancillary_group is not None

    def test_log_file_default(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        tmp_path: Path,
    ):
        """Should set default log file based on output path."""
        from cal_disp.config import InputFileGroup

        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            product_path_group=ProductPathGroup(output_path=tmp_path / "output"),
        )

        workflow = config.to_workflow()

        # Default should be output_path / cal_disp.log
        assert workflow.log_file == tmp_path / "output" / "cal_disp.log"

    def test_custom_log_file(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        tmp_path: Path,
    ):
        """Should accept custom log file."""
        from cal_disp.config import InputFileGroup

        custom_log = tmp_path / "custom.log"

        config = RunConfig(
            input_file_group=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            log_file=custom_log,
        )

        workflow = config.to_workflow()

        assert workflow.log_file == custom_log
