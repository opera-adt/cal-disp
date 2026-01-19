"""Tests for cal-disp config CLI command."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import pytest
import yaml  # type: ignore[import-untyped]
from click.testing import CliRunner

from cal_disp.cli.config import config_cli, create_config


class TestCreateConfig:
    """Tests for create_config function."""

    def test_creates_config_file(
        self,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        sample_frame_id: int,
        sample_grid_version: str,
        sample_grid_type: Literal["constant", "variable"],
    ):
        """Should create a valid configuration file."""
        los_file, dem_file = sample_static_layers
        output_dir = tmp_path / "outputs"
        work_dir = tmp_path / "work"

        config_path = create_config(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=sample_unr_grid_latlon,
            unr_timeseries_dir=sample_unr_timeseries_dir,
            unr_grid_version=sample_grid_version,
            unr_grid_type=sample_grid_type,
            frame_id=sample_frame_id,
            algorithm_params_file=sample_algorithm_params,
            los_file=los_file,
            dem_file=dem_file,
            output_dir=output_dir,
            work_dir=work_dir,
        )

        assert config_path.exists()
        assert config_path.suffix in [".yaml", ".yml"]

    def test_config_contains_required_fields(
        self,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        sample_frame_id: int,
        sample_grid_version: str,
        sample_grid_type: Literal["constant", "variable"],
    ):
        """Generated config should contain all required fields."""
        los_file, dem_file = sample_static_layers
        output_dir = tmp_path / "outputs"
        work_dir = tmp_path / "work"

        config_path = create_config(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=sample_unr_grid_latlon,
            unr_timeseries_dir=sample_unr_timeseries_dir,
            unr_grid_version=sample_grid_version,
            unr_grid_type=sample_grid_type,
            frame_id=sample_frame_id,
            algorithm_params_file=sample_algorithm_params,
            los_file=los_file,
            dem_file=dem_file,
            output_dir=output_dir,
            work_dir=work_dir,
        )

        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "input_file_group" in config["cal_disp_workflow"]
        assert "dynamic_ancillary_group" in config["cal_disp_workflow"]
        assert "product_path_group" in config["cal_disp_workflow"]

    def test_raises_for_missing_required_file(
        self,
        tmp_path: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should raise FileNotFoundError for missing required files."""
        los_file, dem_file = sample_static_layers
        nonexistent_disp = tmp_path / "nonexistent.nc"

        with pytest.raises(FileNotFoundError, match="disp_file"):
            create_config(
                disp_file=nonexistent_disp,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                unr_grid_version="0.2",
                unr_grid_type="constant",
                frame_id=8882,
                algorithm_params_file=sample_algorithm_params,
                los_file=los_file,
                dem_file=dem_file,
                output_dir=tmp_path / "outputs",
                work_dir=tmp_path,
            )

    def test_raises_for_invalid_frame_id(
        self,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should raise ValueError for negative frame ID."""
        los_file, dem_file = sample_static_layers

        with pytest.raises(ValueError, match="frame_id must be >= 0"):
            create_config(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                unr_grid_version="0.2",
                unr_grid_type="constant",
                frame_id=-1,
                algorithm_params_file=sample_algorithm_params,
                los_file=los_file,
                dem_file=dem_file,
                output_dir=tmp_path / "outputs",
                work_dir=tmp_path,
            )

    def test_raises_for_invalid_n_workers(
        self,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should raise ValueError for invalid worker count."""
        los_file, dem_file = sample_static_layers

        with pytest.raises(ValueError, match="n_workers must be >= 1"):
            create_config(
                disp_file=sample_disp_product,
                unr_grid_latlon_file=sample_unr_grid_latlon,
                unr_timeseries_dir=sample_unr_timeseries_dir,
                unr_grid_version="0.2",
                unr_grid_type="constant",
                frame_id=8882,
                algorithm_params_file=sample_algorithm_params,
                los_file=los_file,
                dem_file=dem_file,
                output_dir=tmp_path / "outputs",
                n_workers=0,
                work_dir=tmp_path,
            )

    def test_creates_output_directory(
        self,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should create output directory if it doesn't exist."""
        los_file, dem_file = sample_static_layers
        output_dir = tmp_path / "nonexistent" / "outputs"
        assert not output_dir.exists()

        create_config(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=sample_unr_grid_latlon,
            unr_timeseries_dir=sample_unr_timeseries_dir,
            unr_grid_version="0.2",
            unr_grid_type="constant",
            frame_id=8882,
            algorithm_params_file=sample_algorithm_params,
            los_file=los_file,
            dem_file=dem_file,
            output_dir=output_dir,
            work_dir=tmp_path,
        )

        assert output_dir.exists()

    def test_with_optional_mask_file(
        self,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
        sample_mask_file: Path,
    ):
        """Should accept optional mask file."""
        los_file, dem_file = sample_static_layers

        config_path = create_config(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=sample_unr_grid_latlon,
            unr_timeseries_dir=sample_unr_timeseries_dir,
            unr_grid_version="0.2",
            unr_grid_type="constant",
            frame_id=8882,
            algorithm_params_file=sample_algorithm_params,
            los_file=los_file,
            dem_file=dem_file,
            output_dir=tmp_path / "outputs",
            mask_file=sample_mask_file,
            work_dir=tmp_path,
        )

        with open(config_path) as f:
            config = yaml.safe_load(f)

        assert "mask_file" in config["cal_disp_workflow"]["dynamic_ancillary_group"]

    def test_absolute_paths_by_default(
        self,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should use absolute paths by default."""
        los_file, dem_file = sample_static_layers

        config_path = create_config(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=sample_unr_grid_latlon,
            unr_timeseries_dir=sample_unr_timeseries_dir,
            unr_grid_version="0.2",
            unr_grid_type="constant",
            frame_id=8882,
            algorithm_params_file=sample_algorithm_params,
            los_file=los_file,
            dem_file=dem_file,
            output_dir=tmp_path / "outputs",
            keep_relative=False,
            work_dir=tmp_path,
        )

        with open(config_path) as f:
            config = yaml.safe_load(f)

        disp_path = Path(config["cal_disp_workflow"]["input_file_group"]["disp_file"])
        assert disp_path.is_absolute()


class TestConfigCLI:
    """Tests for config CLI command."""

    def test_basic_invocation(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should run successfully with required arguments."""
        los_file, dem_file = sample_static_layers
        output_dir = tmp_path / "outputs"
        config_file = tmp_path / "runconfig.yaml"

        result = cli_runner.invoke(
            config_cli,
            [
                "--config-file",
                str(config_file),
                "--disp-file",
                str(sample_disp_product),
                "--unr-grid-latlon",
                str(sample_unr_grid_latlon),
                "--unr-grid-dir",
                str(sample_unr_timeseries_dir),
                "--frame-id",
                "8882",
                "--unr-grid-type",
                "constant",
                "--unr-grid-version",
                "0.2",
                "--algorithm-params",
                str(sample_algorithm_params),
                "--los-file",
                str(los_file),
                "--dem-file",
                str(dem_file),
                "--work-dir",
                str(tmp_path),
                "--output-dir",
                str(output_dir),
            ],
        )

        assert result.exit_code == 0
        assert config_file.exists()
        assert "Configuration created" in result.output

    def test_missing_required_option(self, cli_runner: CliRunner):
        """Should fail when required option is missing."""
        result = cli_runner.invoke(config_cli, [])

        assert result.exit_code != 0
        assert "Error" in result.output or "required" in result.output.lower()

    def test_nonexistent_input_file(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should fail gracefully for nonexistent input file."""
        los_file, dem_file = sample_static_layers
        nonexistent = tmp_path / "nonexistent.nc"

        result = cli_runner.invoke(
            config_cli,
            [
                "--disp-file",
                str(nonexistent),
                "--unr-grid-latlon",
                str(sample_unr_grid_latlon),
                "--unr-grid-dir",
                str(sample_unr_timeseries_dir),
                "--frame-id",
                "8882",
                "--unr-grid-type",
                "constant",
                "--unr-grid-version",
                "0.2",
                "--algorithm-params",
                str(sample_algorithm_params),
                "--los-file",
                str(los_file),
                "--dem-file",
                str(dem_file),
                "--output-dir",
                str(tmp_path / "outputs"),
                "--work-dir",
                str(tmp_path),
            ],
        )

        assert result.exit_code != 0

    def test_with_multiple_tropo_files(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should handle multiple tropo files."""
        los_file, dem_file = sample_static_layers

        # Create mock tropo files
        tropo1 = tmp_path / "tropo1.h5"
        tropo2 = tmp_path / "tropo2.h5"
        tropo1.touch()
        tropo2.touch()

        result = cli_runner.invoke(
            config_cli,
            [
                "--disp-file",
                str(sample_disp_product),
                "--unr-grid-latlon",
                str(sample_unr_grid_latlon),
                "--unr-grid-dir",
                str(sample_unr_timeseries_dir),
                "--frame-id",
                "8882",
                "--unr-grid-type",
                "constant",
                "--unr-grid-version",
                "0.2",
                "--algorithm-params",
                str(sample_algorithm_params),
                "--los-file",
                str(los_file),
                "--dem-file",
                str(dem_file),
                "--output-dir",
                str(tmp_path / "outputs"),
                "--ref-tropo-files",
                str(tropo1),
                "--ref-tropo-files",
                str(tropo2),
                "--work-dir",
                str(tmp_path),
            ],
        )

        assert result.exit_code == 0

    def test_custom_worker_settings(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        sample_disp_product: Path,
        sample_unr_grid_latlon: Path,
        sample_unr_timeseries_dir: Path,
        sample_algorithm_params: Path,
        sample_static_layers: tuple[Path, Path],
    ):
        """Should accept custom worker settings."""
        los_file, dem_file = sample_static_layers
        config_file = tmp_path / "runconfig.yaml"

        result = cli_runner.invoke(
            config_cli,
            [
                "--config-file",
                str(config_file),
                "--disp-file",
                str(sample_disp_product),
                "--unr-grid-latlon",
                str(sample_unr_grid_latlon),
                "--unr-grid-dir",
                str(sample_unr_timeseries_dir),
                "--frame-id",
                "8882",
                "--unr-grid-type",
                "constant",
                "--unr-grid-version",
                "0.2",
                "--algorithm-params",
                str(sample_algorithm_params),
                "--los-file",
                str(los_file),
                "--dem-file",
                str(dem_file),
                "--work-dir",
                str(tmp_path),
                "--output-dir",
                str(tmp_path / "outputs"),
                "--n-workers",
                "8",
                "--threads-per-worker",
                "2",
                "--block-shape",
                "256",
                "256",
            ],
        )

        assert result.exit_code == 0

        with open(config_file) as f:
            config = yaml.safe_load(f)

        worker_settings = config["cal_disp_workflow"]["worker_settings"]
        assert worker_settings["n_workers"] == 8
        assert worker_settings["threads_per_worker"] == 2
        assert worker_settings["block_shape"] == [256, 256]
