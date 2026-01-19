"""Tests for cal-disp run and validate CLI commands."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock, patch

import pytest
from click.testing import CliRunner

from cal_disp.cli.run import run_cli, run_main
from cal_disp.cli.validate import validate_cli


class TestRunMain:
    """Tests for run_main function."""

    def test_loads_config_and_runs(self, tmp_path: Path):
        """Should load config file and run calibration."""
        config_file = tmp_path / "runconfig.yaml"
        config_content = """
runconfig:
  name: test_run
  groups:
    input_file_group:
      disp_file: /path/to/disp.nc
      frame_id: 8882
"""
        config_file.write_text(config_content)

        with patch(
            "cal_disp.config.pge_runconfig.RunConfig.from_yaml_file"
        ) as mock_load:
            with patch("cal_disp.main.run") as mock_run:
                mock_runconfig = Mock()
                mock_workflow = Mock()
                mock_runconfig.to_workflow.return_value = mock_workflow
                mock_load.return_value = mock_runconfig
                mock_run.return_value = tmp_path / "output.nc"

                result = run_main(config_file=config_file, debug=False)

                mock_load.assert_called_once_with(config_file)
                mock_run.assert_called_once()
                assert result == tmp_path / "output.nc"

    def test_passes_debug_flag(self, tmp_path: Path):
        """Should pass debug flag to run function."""
        config_file = tmp_path / "runconfig.yaml"
        config_file.write_text("runconfig:\n  name: test\n  groups: {}")

        with patch("cal_disp.config.pge_runconfig.RunConfig.from_yaml_file"):
            with patch("cal_disp.main.run") as mock_run:
                run_main(config_file=config_file, debug=True)

                call_kwargs = mock_run.call_args.kwargs
                assert call_kwargs["debug"] is True


class TestRunCLI:
    """Tests for run CLI command."""

    def test_requires_config_file(self, cli_runner: CliRunner):
        """Should require config file argument."""
        result = cli_runner.invoke(run_cli, [])

        assert result.exit_code != 0

    def test_basic_invocation(self, cli_runner: CliRunner, tmp_path: Path):
        """Should run with valid config file."""
        config_file = tmp_path / "runconfig.yaml"
        config_file.write_text("runconfig:\n  name: test\n  groups: {}")

        with patch("cal_disp.cli.run.run_main") as mock_run:
            mock_run.return_value = tmp_path / "output.nc"

            result = cli_runner.invoke(
                run_cli,
                [str(config_file)],
                obj={"debug": False},
            )

            assert result.exit_code == 0
            mock_run.assert_called_once_with(config_file=config_file, debug=False)

    def test_with_debug_flag(self, cli_runner: CliRunner, tmp_path: Path):
        """Should pass debug flag from context."""
        config_file = tmp_path / "runconfig.yaml"
        config_file.write_text("runconfig:\n  name: test\n  groups: {}")

        with patch("cal_disp.cli.run.run_main") as mock_run:
            mock_run.return_value = tmp_path / "output.nc"

            cli_runner.invoke(
                run_cli,
                [str(config_file)],
                obj={"debug": True},
            )

            call_kwargs = mock_run.call_args.kwargs
            assert call_kwargs["debug"] is True

    def test_nonexistent_config_file(self, cli_runner: CliRunner, tmp_path: Path):
        """Should fail for nonexistent config file."""
        nonexistent = tmp_path / "nonexistent.yaml"

        result = cli_runner.invoke(
            run_cli,
            [str(nonexistent)],
            obj={"debug": False},
        )

        assert result.exit_code != 0


class TestValidateCLI:
    """Tests for validate CLI command."""

    def test_requires_two_files(self, cli_runner: CliRunner):
        """Should require reference and test file arguments."""
        result = cli_runner.invoke(validate_cli, [])

        assert result.exit_code != 0

    def test_basic_invocation(self, cli_runner: CliRunner, tmp_path: Path):
        """Should run with two valid files."""
        ref_file = tmp_path / "reference.nc"
        test_file = tmp_path / "test.nc"
        ref_file.touch()
        test_file.touch()

        with patch("cal_disp.validate.compare_cal_products") as mock_compare:
            mock_compare.return_value = True

            result = cli_runner.invoke(
                validate_cli,
                [str(ref_file), str(test_file)],
                obj={"debug": False},
            )

            assert result.exit_code == 0
            mock_compare.assert_called_once()

    def test_with_custom_tolerance(self, cli_runner: CliRunner, tmp_path: Path):
        """Should accept --tolerance option."""
        ref_file = tmp_path / "reference.nc"
        test_file = tmp_path / "test.nc"
        ref_file.touch()
        test_file.touch()

        with patch("cal_disp.validate.compare_cal_products") as mock_compare:
            mock_compare.return_value = True

            cli_runner.invoke(
                validate_cli,
                [str(ref_file), str(test_file), "--tolerance", "1e-5"],
                obj={"debug": False},
            )

            call_kwargs = mock_compare.call_args.kwargs
            assert call_kwargs["tolerance"] == 1e-5

    def test_with_group_selection(self, cli_runner: CliRunner, tmp_path: Path):
        """Should accept --group option."""
        ref_file = tmp_path / "reference.nc"
        test_file = tmp_path / "test.nc"
        ref_file.touch()
        test_file.touch()

        with patch("cal_disp.validate.compare_cal_products") as mock_compare:
            mock_compare.return_value = True

            cli_runner.invoke(
                validate_cli,
                [str(ref_file), str(test_file), "--group", "main"],
                obj={"debug": False},
            )

            call_kwargs = mock_compare.call_args.kwargs
            assert call_kwargs["group"] == "main"

    def test_fails_on_validation_error(self, cli_runner: CliRunner, tmp_path: Path):
        """Should exit with error when validation fails."""
        ref_file = tmp_path / "reference.nc"
        test_file = tmp_path / "test.nc"
        ref_file.touch()
        test_file.touch()

        with patch("cal_disp.validate.compare_cal_products") as mock_compare:
            mock_compare.return_value = False  # Validation failed

            result = cli_runner.invoke(
                validate_cli,
                [str(ref_file), str(test_file)],
                obj={"debug": False},
            )

            assert result.exit_code != 0

    def test_nonexistent_reference_file(self, cli_runner: CliRunner, tmp_path: Path):
        """Should fail for nonexistent reference file."""
        nonexistent = tmp_path / "nonexistent.nc"
        test_file = tmp_path / "test.nc"
        test_file.touch()

        result = cli_runner.invoke(
            validate_cli,
            [str(nonexistent), str(test_file)],
            obj={"debug": False},
        )

        assert result.exit_code != 0

    def test_with_debug_logging(self, cli_runner: CliRunner, tmp_path: Path):
        """Should enable debug logging when debug flag set."""
        ref_file = tmp_path / "reference.nc"
        test_file = tmp_path / "test.nc"
        ref_file.touch()
        test_file.touch()

        with patch("cal_disp.validate.compare_cal_products") as mock_compare:
            with patch("cal_disp._log.setup_logging") as mock_setup:
                mock_compare.return_value = True

                cli_runner.invoke(
                    validate_cli,
                    [str(ref_file), str(test_file)],
                    obj={"debug": True},
                )

                # Should setup logging with DEBUG level
                mock_setup.assert_called_once()
                call_kwargs = mock_setup.call_args.kwargs
                assert call_kwargs["level"] == "DEBUG"

    @pytest.mark.parametrize(
        "group",
        ["main", "auxiliary", "all"],
    )
    def test_all_group_options(
        self,
        cli_runner: CliRunner,
        tmp_path: Path,
        group: str,
    ):
        """Should accept all valid group options."""
        ref_file = tmp_path / "reference.nc"
        test_file = tmp_path / "test.nc"
        ref_file.touch()
        test_file.touch()

        with patch("cal_disp.validate.compare_cal_products") as mock_compare:
            mock_compare.return_value = True

            result = cli_runner.invoke(
                validate_cli,
                [str(ref_file), str(test_file), "--group", group],
                obj={"debug": False},
            )

            assert result.exit_code == 0

    def test_invalid_group_option(self, cli_runner: CliRunner, tmp_path: Path):
        """Should reject invalid group option."""
        ref_file = tmp_path / "reference.nc"
        test_file = tmp_path / "test.nc"
        ref_file.touch()
        test_file.touch()

        result = cli_runner.invoke(
            validate_cli,
            [str(ref_file), str(test_file), "--group", "invalid"],
            obj={"debug": False},
        )

        assert result.exit_code != 0
