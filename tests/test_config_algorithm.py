"""Tests for algorithm parameters configuration."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml  # type: ignore[import-untyped]

from cal_disp.config._algorithm import (
    AlgorithmParameters,
    CalibrationMethod,
    CalibrationOptions,
    SavitskyGoleyOptions,
    UnwrapCorrection,
)


class TestCalibrationMethod:
    """Tests for CalibrationMethod enum."""

    def test_has_expected_values(self):
        """Should have all expected calibration methods."""
        assert CalibrationMethod.GAUSSIAN_FFT == "gaussian_fft"
        assert CalibrationMethod.HANNING_FFT == "hanning_fft"
        assert CalibrationMethod.SAVITSKY_GOLEY == "savitsky_goley"

    def test_enum_membership(self):
        """Should validate enum membership."""
        assert "gaussian_fft" in [m.value for m in CalibrationMethod]
        assert "invalid" not in [m.value for m in CalibrationMethod]


class TestUnwrapCorrection:
    """Tests for UnwrapCorrection configuration."""

    def test_default_value(self):
        """Should default to True."""
        config = UnwrapCorrection()
        assert config.run_unwrap_correction is True

    def test_explicit_false(self):
        """Should accept explicit False."""
        config = UnwrapCorrection(run_unwrap_correction=False)
        assert config.run_unwrap_correction is False

    def test_yaml_serialization(self, tmp_path: Path):
        """Should serialize to YAML."""
        config = UnwrapCorrection(run_unwrap_correction=False)
        yaml_file = tmp_path / "unwrap.yaml"

        config.to_yaml(yaml_file)

        with open(yaml_file) as f:
            data = yaml.safe_load(f)
        assert data["run_unwrap_correction"] is False

    def test_yaml_deserialization(self, tmp_path: Path):
        """Should deserialize from YAML."""
        yaml_file = tmp_path / "unwrap.yaml"
        yaml_file.write_text("run_unwrap_correction: false\n")

        config = UnwrapCorrection.from_yaml(yaml_file)
        assert config.run_unwrap_correction is False


class TestSavitskyGoleyOptions:
    """Tests for Savitsky-Golay filter options."""

    def test_all_none_by_default(self):
        """All fields should default to None."""
        options = SavitskyGoleyOptions()

        assert options.window_x_size is None
        assert options.window_y_size is None
        assert options.window_overlap_x_size is None
        assert options.window_overlap_y_size is None
        assert options.window_extend_x_size is None
        assert options.window_extend_y_size is None

    def test_valid_window_sizes(self):
        """Should accept positive window sizes."""
        options = SavitskyGoleyOptions(
            window_x_size=100,
            window_y_size=100,
        )

        assert options.window_x_size == 100
        assert options.window_y_size == 100

    def test_rejects_zero_window_size(self):
        """Should reject zero window sizes."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            SavitskyGoleyOptions(window_x_size=0)

    def test_rejects_negative_window_size(self):
        """Should reject negative window sizes."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            SavitskyGoleyOptions(window_y_size=-1)

    def test_rejects_negative_overlap(self):
        """Should reject negative overlap."""
        with pytest.raises(ValueError, match="greater than or equal to 0"):
            SavitskyGoleyOptions(window_overlap_x_size=-1)

    def test_rejects_negative_extend(self):
        """Should reject negative extend."""
        with pytest.raises(ValueError, match="greater than or equal to 0"):
            SavitskyGoleyOptions(window_extend_y_size=-5)

    def test_all_fields_specified(self):
        """Should accept all fields when specified."""
        options = SavitskyGoleyOptions(
            window_x_size=100,
            window_y_size=100,
            window_overlap_x_size=50,
            window_overlap_y_size=50,
            window_extend_x_size=10,
            window_extend_y_size=10,
        )

        assert options.window_x_size == 100
        assert options.window_overlap_x_size == 50
        assert options.window_extend_x_size == 10


class TestCalibrationOptions:
    """Tests for CalibrationOptions configuration."""

    def test_defaults(self):
        """Should have expected default values."""
        options = CalibrationOptions()

        assert options.cal_method == CalibrationMethod.SAVITSKY_GOLEY
        assert options.run_interpolation is True
        assert options.run_downsample is True
        assert options.downsample_factor == 10

    def test_custom_method(self):
        """Should accept different calibration methods."""
        options = CalibrationOptions(cal_method=CalibrationMethod.GAUSSIAN_FFT)

        assert options.cal_method == CalibrationMethod.GAUSSIAN_FFT

    def test_disable_interpolation(self):
        """Should allow disabling interpolation."""
        options = CalibrationOptions(run_interpolation=False)

        assert options.run_interpolation is False

    def test_custom_downsample_factor(self):
        """Should accept custom downsample factor."""
        options = CalibrationOptions(downsample_factor=5)

        assert options.downsample_factor == 5

    def test_rejects_zero_downsample(self):
        """Should reject zero downsample factor."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            CalibrationOptions(downsample_factor=0)

    def test_rejects_negative_downsample(self):
        """Should reject negative downsample factor."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            CalibrationOptions(downsample_factor=-1)


class TestAlgorithmParameters:
    """Tests for complete AlgorithmParameters."""

    def test_defaults(self):
        """Should create with all default values."""
        params = AlgorithmParameters()

        assert params.unwrap_correction.run_unwrap_correction is True
        assert params.calibration_options.cal_method == CalibrationMethod.SAVITSKY_GOLEY
        assert params.savitsky_goley_options.window_x_size is None

    def test_create_default_classmethod(self):
        """Should create defaults via classmethod."""
        params = AlgorithmParameters.create_default()

        assert isinstance(params, AlgorithmParameters)
        assert params.calibration_options.run_interpolation is True

    def test_create_example_classmethod(self):
        """Should create example with all fields specified."""
        params = AlgorithmParameters.create_example()

        assert params.unwrap_correction.run_unwrap_correction is True
        assert params.calibration_options.downsample_factor == 10
        assert params.savitsky_goley_options.window_x_size == 100
        assert params.savitsky_goley_options.window_y_size == 100
        assert params.savitsky_goley_options.window_overlap_x_size == 50

    def test_custom_configuration(self):
        """Should accept custom nested configuration."""
        params = AlgorithmParameters(
            unwrap_correction=UnwrapCorrection(run_unwrap_correction=False),
            calibration_options=CalibrationOptions(
                cal_method=CalibrationMethod.HANNING_FFT,
                downsample_factor=5,
            ),
            savitsky_goley_options=SavitskyGoleyOptions(
                window_x_size=200,
                window_y_size=200,
            ),
        )

        assert params.unwrap_correction.run_unwrap_correction is False
        assert params.calibration_options.cal_method == CalibrationMethod.HANNING_FFT
        assert params.savitsky_goley_options.window_x_size == 200

    def test_yaml_round_trip(self, tmp_path: Path):
        """Should serialize and deserialize via YAML."""
        params = AlgorithmParameters.create_example()
        yaml_file = tmp_path / "algorithm_params.yaml"

        # Save
        params.to_yaml(yaml_file)
        assert yaml_file.exists()

        # Load
        loaded = AlgorithmParameters.from_yaml(yaml_file)

        assert (
            loaded.calibration_options.cal_method
            == params.calibration_options.cal_method
        )
        assert (
            loaded.savitsky_goley_options.window_x_size
            == params.savitsky_goley_options.window_x_size
        )

    def test_partial_savitsky_options(self):
        """Should allow partial Savitsky-Golay options."""
        params = AlgorithmParameters(
            savitsky_goley_options=SavitskyGoleyOptions(
                window_x_size=100,
                # Leave other fields as None
            )
        )

        assert params.savitsky_goley_options.window_x_size == 100
        assert params.savitsky_goley_options.window_y_size is None
        assert params.savitsky_goley_options.window_overlap_x_size is None

    def test_nested_validation_propagates(self):
        """Should propagate validation errors from nested models."""
        with pytest.raises(ValueError):
            AlgorithmParameters(
                calibration_options=CalibrationOptions(downsample_factor=-1)
            )

    def test_forbids_extra_fields(self):
        """Should forbid extra fields not in schema."""
        with pytest.raises(ValueError, match="extra"):
            AlgorithmParameters(invalid_field="test")

    def test_yaml_with_comments(self, tmp_path: Path):
        """Should generate YAML with field descriptions."""
        params = AlgorithmParameters.create_default()
        yaml_file = tmp_path / "params_commented.yaml"

        params.to_yaml(yaml_file, with_comments=True)

        content = yaml_file.read_text()
        # Comments should exist
        assert "#" in content

    def test_modifying_after_creation(self):
        """Should allow modifying fields after creation."""
        params = AlgorithmParameters()

        params.calibration_options.downsample_factor = 20

        assert params.calibration_options.downsample_factor == 20

    @pytest.mark.parametrize(
        "method",
        [
            CalibrationMethod.GAUSSIAN_FFT,
            CalibrationMethod.HANNING_FFT,
            CalibrationMethod.SAVITSKY_GOLEY,
        ],
    )
    def test_all_calibration_methods(self, method: CalibrationMethod):
        """Should work with all calibration methods."""
        params = AlgorithmParameters(
            calibration_options=CalibrationOptions(cal_method=method)
        )

        assert params.calibration_options.cal_method == method
