"""Every PGE runconfig field is either applied or rejected/warned about."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from unittest.mock import patch

import pytest
from click.testing import CliRunner
from pydantic import ValidationError

from cal_disp.cli.config import config_cli, create_config
from cal_disp.config import (
    DynamicAncillaryFileGroup,
    InputFileGroup,
    StaticAncillaryFileGroup,
)
from cal_disp.config._algorithm import AlgorithmParameters
from cal_disp.config.pge_runconfig import (
    OutputOptions,
    PrimaryExecutable,
    ProductPathGroup,
    RunConfig,
)
from cal_disp.config.workflow import CalibrationWorkflow


@pytest.fixture
def groups(
    sample_disp_product: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
) -> dict:
    """Keyword arguments for the two required file groups."""
    lookup_file, tenv8_dir = sample_unr_data
    return {
        "input": {
            "disp_file": sample_disp_product,
            "frame_id": 8882,
            "unr_grid_latlon_file": lookup_file,
            "unr_timeseries_dir": tenv8_dir,
            "unr_grid_version": "0.2",
            "unr_grid_type": "constant",
        },
        "dynamic": {
            "algorithm_parameters_file": sample_algorithm_params,
            "static_los_file": sample_static_los,
            "static_dem_file": sample_static_dem,
        },
    }


def _workflow(
    tmp_path: Path, groups: dict, static=None, **dynamic
) -> CalibrationWorkflow:
    return CalibrationWorkflow(
        input_options=InputFileGroup(**groups["input"]),
        dynamic_ancillary_options=DynamicAncillaryFileGroup(
            **groups["dynamic"], **dynamic
        ),
        static_ancillary_options=static,
        work_directory=tmp_path / "work",
        output_directory=tmp_path / "out",
    )


def _run_main(workflow: CalibrationWorkflow, tmp_path: Path):
    """Run ``main.run`` with the calibration itself replaced by a recorder."""
    from cal_disp import main

    with (
        patch.object(main, "run_calibration") as run_calibration,
        patch.object(main, "make_browse_image_from_nc"),
    ):
        run_calibration.return_value = tmp_path / "out" / "product.nc"
        main.run(workflow)
    return run_calibration.call_args.kwargs


# algorithm_parameters_overrides_json is applied


class TestAlgorithmOverrides:
    def test_flat_keys_override_calibration_options(self):
        params = AlgorithmParameters()
        new = params.with_overrides({"downsample_factor": 3, "grid_type": "variable"})

        assert new.calibration_options.downsample_factor == 3
        assert new.calibration_options.grid_type == "variable"
        # Everything else, and the original, is unchanged
        assert new.calibration_options.window_size_meters == (
            params.calibration_options.window_size_meters
        )
        assert params.calibration_options.downsample_factor != 3

    def test_nested_groups_are_merged(self):
        params = AlgorithmParameters()
        new = params.with_overrides(
            {
                "calibration_options": {
                    "posting_meters": 90.0,
                    "savitzky_golay": {"window_length": 31},
                }
            }
        )

        assert new.calibration_options.posting_meters == 90.0
        assert new.calibration_options.savitzky_golay.window_length == 31
        assert new.calibration_options.savitzky_golay.polyorder == (
            params.calibration_options.savitzky_golay.polyorder
        )

    def test_empty_overrides_return_self(self):
        params = AlgorithmParameters()
        assert params.with_overrides({}) is params

    def test_unknown_option_raises(self):
        with pytest.raises(ValidationError, match="cal_method"):
            AlgorithmParameters().with_overrides({"cal_method": "hanning_fft"})

    def test_invalid_value_raises(self):
        with pytest.raises(ValidationError, match="grid_type"):
            AlgorithmParameters().with_overrides({"grid_type": "sometimes"})

    def test_main_run_applies_overrides_for_the_frame(self, tmp_path, groups, caplog):
        overrides_file = tmp_path / "overrides.json"
        overrides_file.write_text(
            json.dumps(
                {
                    "data": {
                        "8882": {"downsample_factor": 2, "posting_meters": 50.0},
                        "1234": {"downsample_factor": 7},
                    }
                }
            )
        )
        static = StaticAncillaryFileGroup(
            algorithm_parameters_overrides_json=overrides_file
        )
        with caplog.at_level(logging.INFO, logger="cal_disp"):
            kwargs = _run_main(_workflow(tmp_path, groups, static=static), tmp_path)

        cal = kwargs["algorithm_parameters"].calibration_options
        assert cal.downsample_factor == 2  # the algorithm file has 1
        assert cal.posting_meters == 50.0  # the algorithm file has 100
        assert cal.window_size_meters == 5000.0  # untouched, from the file
        assert "Algorithm parameter overrides for frame 8882" in caplog.text

    def test_main_run_without_entry_for_the_frame(self, tmp_path, groups, caplog):
        overrides_file = tmp_path / "overrides.json"
        overrides_file.write_text(json.dumps({"1234": {"downsample_factor": 7}}))
        static = StaticAncillaryFileGroup(
            algorithm_parameters_overrides_json=overrides_file
        )
        with caplog.at_level(logging.INFO, logger="cal_disp"):
            kwargs = _run_main(_workflow(tmp_path, groups, static=static), tmp_path)

        assert kwargs["algorithm_parameters"].calibration_options.downsample_factor == 1
        assert "No algorithm parameter overrides for frame 8882" in caplog.text

    def test_main_run_rejects_unknown_override(self, tmp_path, groups):
        overrides_file = tmp_path / "overrides.json"
        overrides_file.write_text(json.dumps({"8882": {"no_such_option": 1}}))
        static = StaticAncillaryFileGroup(
            algorithm_parameters_overrides_json=overrides_file
        )
        with pytest.raises(ValidationError, match="no_such_option"):
            _run_main(_workflow(tmp_path, groups, static=static), tmp_path)


# Fields that are not supported in this release fail loudly


class TestUnsupportedFields:
    @pytest.mark.parametrize("field", ["iono_files", "tiles_files"])
    def test_non_empty_file_list_raises(self, groups, tmp_path, field):
        extra = tmp_path / "extra.nc"
        extra.touch()
        with pytest.raises(ValidationError, match=f"{field} is not supported"):
            DynamicAncillaryFileGroup(**groups["dynamic"], **{field: [extra]})

    @pytest.mark.parametrize("field", ["iono_files", "tiles_files"])
    @pytest.mark.parametrize("value", [None, []])
    def test_empty_file_list_is_accepted(self, groups, field, value):
        group = DynamicAncillaryFileGroup(**groups["dynamic"], **{field: value})
        assert getattr(group, field) == []

    def test_assigning_a_file_list_raises(self, groups, tmp_path):
        group = DynamicAncillaryFileGroup(**groups["dynamic"])
        with pytest.raises(ValidationError, match="iono_files is not supported"):
            group.iono_files = [tmp_path / "iono.nc"]

    def test_output_format_other_than_netcdf_raises(self):
        assert OutputOptions(output_format="netcdf").output_format == "netcdf"
        with pytest.raises(ValidationError, match="output_format='hdf5' is not"):
            OutputOptions(output_format="hdf5")

    def test_runconfig_yaml_with_unsupported_field_fails_to_load(
        self, groups, tmp_path
    ):
        config = RunConfig(
            input_file_group=InputFileGroup(**groups["input"]),
            dynamic_ancillary_group=DynamicAncillaryFileGroup(**groups["dynamic"]),
        )
        yaml_file = tmp_path / "runconfig.yaml"
        config.to_yaml(yaml_file)
        assert (
            RunConfig.from_yaml_file(yaml_file).output_options.output_format == "netcdf"
        )

        text = yaml_file.read_text()
        assert "output_format: netcdf" in text
        yaml_file.write_text(
            text.replace("output_format: netcdf", "output_format: hdf5")
        )
        with pytest.raises(ValidationError, match="not supported in this release"):
            RunConfig.from_yaml_file(yaml_file)

    def test_create_config_rejects_iono_files(self, groups, tmp_path):
        iono = tmp_path / "iono.nc"
        iono.touch()
        with pytest.raises(ValueError, match="iono_files is not supported"):
            create_config(
                disp_file=groups["input"]["disp_file"],
                frame_id=8882,
                unr_grid_latlon_file=groups["input"]["unr_grid_latlon_file"],
                unr_timeseries_dir=groups["input"]["unr_timeseries_dir"],
                unr_grid_version="0.2",
                unr_grid_type="constant",
                algorithm_params_file=groups["dynamic"]["algorithm_parameters_file"],
                los_file=groups["dynamic"]["static_los_file"],
                dem_file=groups["dynamic"]["static_dem_file"],
                output_dir=tmp_path / "out",
                work_dir=tmp_path / "work",
                iono_files=[iono],
            )


# Fields that are accepted but unused are reported with a WARNING


class TestIgnoredFieldWarnings:
    def test_mask_file_warns(self, tmp_path, groups, sample_mask_file, caplog):
        workflow = _workflow(tmp_path, groups, mask_file=sample_mask_file)
        with caplog.at_level(logging.WARNING, logger="cal_disp"):
            _run_main(workflow, tmp_path)

        assert "mask_file" in caplog.text
        assert "is not applied in this release" in caplog.text

    def test_no_mask_file_no_warning(self, tmp_path, groups, caplog):
        with caplog.at_level(logging.WARNING, logger="cal_disp"):
            _run_main(_workflow(tmp_path, groups), tmp_path)

        assert "mask_file" not in caplog.text

    def _runconfig(self, groups, **kwargs) -> RunConfig:
        return RunConfig(
            input_file_group=InputFileGroup(**groups["input"]),
            dynamic_ancillary_group=DynamicAncillaryFileGroup(**groups["dynamic"]),
            **kwargs,
        )

    def test_product_path_differing_from_output_path_warns(
        self, tmp_path, groups, caplog
    ):
        config = self._runconfig(
            groups,
            product_path_group=ProductPathGroup(
                product_path=tmp_path / "pge_products",
                scratch_path=tmp_path / "scratch",
                sas_output_path=tmp_path / "output",
            ),
        )
        with caplog.at_level(logging.WARNING, logger="cal_disp"):
            workflow = config.to_workflow()

        assert "product_path_group.product_path" in caplog.text
        assert "is not used by the SAS" in caplog.text
        assert workflow.output_directory == (tmp_path / "output").resolve()

    def test_custom_product_type_warns(self, tmp_path, groups, caplog):
        config = self._runconfig(
            groups,
            primary_executable=PrimaryExecutable(product_type="CUSTOM"),
            product_path_group=ProductPathGroup(
                product_path=tmp_path / "output",
                scratch_path=tmp_path / "scratch",
                sas_output_path=tmp_path / "output",
            ),
        )
        with caplog.at_level(logging.WARNING, logger="cal_disp"):
            config.to_workflow()

        assert "primary_executable.product_type='CUSTOM' is not used" in caplog.text
        assert "product_path" not in caplog.text

    def test_default_runconfig_does_not_warn(self, tmp_path, groups, caplog):
        config = self._runconfig(
            groups,
            product_path_group=ProductPathGroup(
                product_path=tmp_path / "output",
                scratch_path=tmp_path / "scratch",
                sas_output_path=tmp_path / "output",
            ),
        )
        with caplog.at_level(logging.WARNING, logger="cal_disp"):
            config.to_workflow()

        assert caplog.text == ""


# `cal-disp config -c PATH` honours the full path


class TestConfigFilePath:
    def _args(self, groups, tmp_path) -> list[str]:
        return [
            "--disp-file",
            str(groups["input"]["disp_file"]),
            "--unr-grid-latlon",
            str(groups["input"]["unr_grid_latlon_file"]),
            "--unr-grid-dir",
            str(groups["input"]["unr_timeseries_dir"]),
            "--frame-id",
            "8882",
            "--unr-grid-version",
            "0.2",
            "--algorithm-params",
            str(groups["dynamic"]["algorithm_parameters_file"]),
            "--los-file",
            str(groups["dynamic"]["static_los_file"]),
            "--dem-file",
            str(groups["dynamic"]["static_dem_file"]),
            "--output-dir",
            str(tmp_path / "out"),
            "--work-dir",
            str(tmp_path / "work"),
        ]

    def test_config_file_directory_is_honoured(self, groups, tmp_path):
        config_file = tmp_path / "configs" / "nested" / "my_runconfig.yaml"
        result = CliRunner().invoke(
            config_cli, [*self._args(groups, tmp_path), "-c", str(config_file)]
        )

        assert result.exit_code == 0, result.output
        assert config_file.is_file()
        assert f"Configuration created: {config_file}" in result.output
        # Not written into the work directory under the same name
        assert not (tmp_path / "work" / "my_runconfig.yaml").exists()
        assert RunConfig.from_yaml_file(config_file).input_file_group.frame_id == 8882

    def test_default_is_runconfig_yaml_in_the_work_dir(self, groups, tmp_path):
        result = CliRunner().invoke(config_cli, self._args(groups, tmp_path))

        assert result.exit_code == 0, result.output
        assert (tmp_path / "work" / "runconfig.yaml").is_file()

    def test_relative_config_file_is_relative_to_cwd(
        self, groups, tmp_path, monkeypatch
    ):
        cwd = tmp_path / "cwd"
        cwd.mkdir()
        monkeypatch.chdir(cwd)
        result = CliRunner().invoke(
            config_cli, [*self._args(groups, tmp_path), "-c", "configs/rc.yaml"]
        )

        assert result.exit_code == 0, result.output
        assert (cwd / "configs" / "rc.yaml").is_file()
        assert not (tmp_path / "work" / "rc.yaml").exists()

    @pytest.mark.parametrize(
        "extra",
        [
            ["--iono-files", "somefile"],
            ["--tiles-files", "somefile"],
            ["--output-format", "hdf5"],
        ],
    )
    def test_unsupported_cli_options_are_rejected(self, groups, tmp_path, extra):
        result = CliRunner().invoke(config_cli, [*self._args(groups, tmp_path), *extra])

        assert result.exit_code == 2
        assert not (tmp_path / "work" / "runconfig.yaml").exists()
