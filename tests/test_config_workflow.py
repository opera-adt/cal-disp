"""Tests for CalibrationWorkflow configuration."""

from __future__ import annotations

from pathlib import Path

from cal_disp.config.workflow import CalibrationWorkflow


class TestCalibrationWorkflowBasics:
    """Tests for basic CalibrationWorkflow functionality."""

    def test_creates_with_minimal_config(self):
        """Should create with minimal configuration."""
        workflow = CalibrationWorkflow()

        assert workflow.work_directory.is_absolute()
        assert workflow.output_directory.is_absolute()
        assert workflow.input_options is None

    def test_creates_with_relative_paths(self):
        """Should keep paths relative when requested."""
        workflow = CalibrationWorkflow(keep_paths_relative=True)
        assert workflow.work_directory == Path()
        assert workflow.output_directory == Path()
        assert workflow.input_options is None

    def test_creates_with_custom_directories(self, tmp_path: Path):
        """Should accept custom directories."""
        work_dir = tmp_path / "work"
        output_dir = tmp_path / "output"

        workflow = CalibrationWorkflow(
            work_directory=work_dir,
            output_directory=output_dir,
        )

        assert workflow.work_directory == work_dir
        assert workflow.output_directory == output_dir

    def test_default_log_file_location(self, tmp_path: Path):
        """Should set default log file in work directory."""
        work_dir = tmp_path / "work"

        workflow = CalibrationWorkflow(work_directory=work_dir)

        assert workflow.log_file == work_dir / "cal_disp.log"

    def test_custom_log_file(self, tmp_path: Path):
        """Should accept custom log file."""
        custom_log = tmp_path / "custom.log"

        workflow = CalibrationWorkflow(log_file=custom_log)

        assert workflow.log_file == custom_log

    def test_keep_paths_relative_default(self):
        """Should resolve to absolute paths by default."""
        workflow = CalibrationWorkflow(
            work_directory=Path("./work"),
            output_directory=Path("./output"),
        )

        # Paths should be absolute
        assert workflow.work_directory.is_absolute()
        assert workflow.output_directory.is_absolute()

    def test_keep_paths_relative_true(self):
        """Should keep paths relative when requested."""
        workflow = CalibrationWorkflow(
            work_directory=Path("./work"),
            output_directory=Path("./output"),
            keep_paths_relative=True,
        )

        # Paths should remain relative
        assert not workflow.work_directory.is_absolute()
        assert not workflow.output_directory.is_absolute()


class TestCalibrationWorkflowWithInputs:
    """Tests for workflow with input configurations."""

    def test_with_input_options(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should accept input options."""
        from cal_disp.config import InputFileGroup

        input_opts = InputFileGroup(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=sample_unr_grid_latlon,
            unr_timeseries_dir=sample_unr_timeseries_dir,
            frame_id=8882,
            unr_grid_version="0.2",
            unr_grid_type="constant",
        )

        workflow = CalibrationWorkflow(input_options=input_opts)

        assert workflow.input_options is not None
        assert workflow.input_options.frame_id == 8882

    def test_with_dynamic_ancillary(
        self,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should accept dynamic ancillary files."""
        from cal_disp.config import DynamicAncillaryFileGroup

        los_file, dem_file = sample_static_layers

        dynamic_opts = DynamicAncillaryFileGroup(
            algorithm_parameters_file=sample_algorithm_params,
            static_los_file=los_file,
            static_dem_file=dem_file,
        )

        workflow = CalibrationWorkflow(dynamic_ancillary_options=dynamic_opts)

        assert workflow.dynamic_ancillary_options is not None

    def test_with_static_ancillary(self, tmp_path: Path):
        """Should accept static ancillary files."""
        from cal_disp.config import StaticAncillaryFileGroup

        override_file = tmp_path / "overrides.json"
        override_file.write_text("{}")

        static_opts = StaticAncillaryFileGroup(
            algorithm_parameters_overrides_json=override_file
        )

        workflow = CalibrationWorkflow(static_ancillary_options=static_opts)

        assert workflow.static_ancillary_options is not None
        assert workflow.static_ancillary_options.has_algorithm_overrides()


class TestCalibrationWorkflowValidation:
    """Tests for workflow validation."""

    def test_validate_ready_without_inputs(self):
        """Should report not ready without inputs."""
        workflow = CalibrationWorkflow()

        validation = workflow.validate_ready_to_run()

        assert not validation["ready"]
        assert len(validation["errors"]) > 0

    def test_validate_ready_missing_dynamic_ancillary(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
    ):
        """Should report error for missing dynamic ancillary."""
        from cal_disp.config import InputFileGroup

        workflow = CalibrationWorkflow(
            input_options=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )
        )

        validation = workflow.validate_ready_to_run()

        assert not validation["ready"]
        assert any("dynamic_ancillary" in e for e in validation["errors"])

    def test_validate_ready_complete(
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

        workflow = CalibrationWorkflow(
            input_options=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            dynamic_ancillary_options=DynamicAncillaryFileGroup(
                algorithm_parameters_file=sample_algorithm_params,
                static_los_file=los_file,
                static_dem_file=dem_file,
            ),
        )

        validation = workflow.validate_ready_to_run()

        assert validation["ready"]
        assert len(validation["errors"]) == 0

    def test_validate_input_files_exist(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should validate file existence."""
        from cal_disp.config import DynamicAncillaryFileGroup, InputFileGroup

        los_file, dem_file = sample_static_layers

        workflow = CalibrationWorkflow(
            input_options=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            dynamic_ancillary_options=DynamicAncillaryFileGroup(
                algorithm_parameters_file=sample_algorithm_params,
                static_los_file=los_file,
                static_dem_file=dem_file,
            ),
        )

        results = workflow.validate_input_files_exist()

        # All files should exist
        assert all(info["exists"] for info in results.values())

    def test_get_missing_files_none(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should return empty list when all files exist."""
        from cal_disp.config import DynamicAncillaryFileGroup, InputFileGroup

        los_file, dem_file = sample_static_layers

        workflow = CalibrationWorkflow(
            input_options=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            dynamic_ancillary_options=DynamicAncillaryFileGroup(
                algorithm_parameters_file=sample_algorithm_params,
                static_los_file=los_file,
                static_dem_file=dem_file,
            ),
        )

        missing = workflow.get_missing_files()

        assert len(missing) == 0

    def test_get_missing_files_some(self, tmp_path: Path):
        """Should report missing files."""
        from cal_disp.config import InputFileGroup

        nonexistent = tmp_path / "nonexistent.nc"

        # create a dummy .tenv8 file
        (tmp_path / "dummy.tenv8").write_text("fake content")

        workflow = CalibrationWorkflow(
            input_options=InputFileGroup(
                disp_file=nonexistent,
                unr_grid_latlon_file=tmp_path / "grid_latlon_lookup_v0.2.txt",
                unr_timeseries_dir=tmp_path,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            )
        )

        missing = workflow.get_missing_files()

        assert len(missing) > 0


class TestCalibrationWorkflowDirectories:
    """Tests for directory creation and management."""

    def test_create_directories(self, tmp_path: Path):
        """Should create work and output directories."""
        work_dir = tmp_path / "work"
        output_dir = tmp_path / "output"

        workflow = CalibrationWorkflow(
            work_directory=work_dir,
            output_directory=output_dir,
        )

        assert not work_dir.exists()
        assert not output_dir.exists()

        workflow.create_directories()

        assert work_dir.exists()
        assert output_dir.exists()

    def test_create_directories_with_log(self, tmp_path: Path):
        """Should create parent directory for log file."""
        log_dir = tmp_path / "logs"
        log_file = log_dir / "custom.log"

        workflow = CalibrationWorkflow(log_file=log_file)

        assert not log_dir.exists()

        workflow.create_directories()

        assert log_dir.exists()

    def test_create_directories_exist_ok(self, tmp_path: Path):
        """Should not error if directories already exist."""
        work_dir = tmp_path / "work"
        work_dir.mkdir()

        workflow = CalibrationWorkflow(work_directory=work_dir)

        # Should not raise
        workflow.create_directories(exist_ok=True)


class TestCalibrationWorkflowLogging:
    """Tests for logging setup."""

    def test_setup_logging_default(self, tmp_path: Path):
        """Should setup logging with defaults."""
        workflow = CalibrationWorkflow(
            work_directory=tmp_path / "work",
            log_file=tmp_path / "test.log",
        )

        logger = workflow.setup_logging()

        assert logger.name == "cal_disp"
        assert len(logger.handlers) == 2  # Console + file

    def test_setup_logging_custom_level(self, tmp_path: Path):
        """Should accept custom log level."""
        import logging

        workflow = CalibrationWorkflow(
            work_directory=tmp_path / "work",
            log_file=tmp_path / "test.log",
        )

        logger = workflow.setup_logging(level=logging.DEBUG)

        assert logger.level == logging.DEBUG

    def test_setup_logging_creates_log_dir(self, tmp_path: Path):
        """Should create log file directory if needed."""
        log_dir = tmp_path / "logs"
        log_file = log_dir / "test.log"

        workflow = CalibrationWorkflow(log_file=log_file)

        assert not log_dir.exists()

        workflow.setup_logging()

        assert log_dir.exists()


class TestCalibrationWorkflowSummary:
    """Tests for summary generation."""

    def test_summary_minimal(self):
        """Should generate summary for minimal config."""
        workflow = CalibrationWorkflow()

        summary = workflow.summary()

        assert "Calibration Workflow Configuration" in summary
        assert "NOT READY" in summary or "Workflow Status" in summary

    def test_summary_complete(
        self,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should generate detailed summary for complete config."""
        from cal_disp.config import DynamicAncillaryFileGroup, InputFileGroup

        los_file, dem_file = sample_static_layers

        workflow = CalibrationWorkflow(
            input_options=InputFileGroup(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                frame_id=8882,
                unr_grid_version="0.2",
                unr_grid_type="constant",
            ),
            dynamic_ancillary_options=DynamicAncillaryFileGroup(
                algorithm_parameters_file=sample_algorithm_params,
                static_los_file=los_file,
                static_dem_file=dem_file,
            ),
        )

        summary = workflow.summary()

        assert "8882" in summary  # Frame ID
        assert "READY" in summary
        assert "Worker Settings" in summary


class TestCalibrationWorkflowClassMethods:
    """Tests for factory class methods."""

    def test_create_example(self):
        """Should create example configuration."""
        workflow = CalibrationWorkflow.create_example()

        assert workflow.input_options is not None
        assert workflow.input_options.frame_id == 8882
        assert workflow.dynamic_ancillary_options is not None

    def test_create_minimal(self):
        """Should create minimal configuration."""
        workflow = CalibrationWorkflow.create_minimal()

        assert workflow.input_options is None
        assert workflow.work_directory == Path("./work")
        assert workflow.output_directory == Path("./output")

    def test_example_keeps_paths_relative(self):
        """Example should use relative paths."""
        workflow = CalibrationWorkflow.create_example()

        # Example has keep_paths_relative=True
        assert workflow.keep_paths_relative is True
