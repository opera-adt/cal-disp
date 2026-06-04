"""Tests for algorithm parameters configuration."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml  # type: ignore[import-untyped]

from cal_disp.config._algorithm import (
    AlgorithmParameters,
    CalibrationOptions,
    FFTFilterOptions,
    SavitzkyGolayOptions,
)


class TestSavitzkyGolayOptions:
    """Tests for Savitzky-Golay filter options."""

    def test_defaults(self):
        """Should have expected default values."""
        options = SavitzkyGolayOptions()

        assert options.window_length == 51
        assert options.polyorder == 3
        assert options.deriv == 0

    def test_valid_window_length(self):
        """Should accept valid window lengths."""
        options = SavitzkyGolayOptions(window_length=101)

        assert options.window_length == 101

    def test_rejects_window_length_below_3(self):
        """Should reject window lengths below 3."""
        with pytest.raises(ValueError, match="greater than or equal to 3"):
            SavitzkyGolayOptions(window_length=2)

    def test_rejects_negative_polyorder(self):
        """Should reject negative polynomial order."""
        with pytest.raises(ValueError, match="greater than or equal to 0"):
            SavitzkyGolayOptions(polyorder=-1)

    def test_rejects_negative_deriv(self):
        """Should reject negative derivative order."""
        with pytest.raises(ValueError, match="greater than or equal to 0"):
            SavitzkyGolayOptions(deriv=-1)

    def test_all_fields_specified(self):
        """Should accept all fields when specified."""
        options = SavitzkyGolayOptions(window_length=99, polyorder=5, deriv=1)

        assert options.window_length == 99
        assert options.polyorder == 5
        assert options.deriv == 1


class TestFFTFilterOptions:
    """Tests for FFT filter options."""

    def test_defaults(self):
        """Should have expected default values."""
        options = FFTFilterOptions()

        assert options.gaussian_sigma == pytest.approx(0.1)
        assert options.butterworth_order == 4
        assert options.spatial_domain is False
        assert options.taper_edges is True
        assert options.taper_width == pytest.approx(0.05)

    def test_rejects_non_positive_sigma(self):
        """Should reject non-positive gaussian_sigma."""
        with pytest.raises(ValueError):
            FFTFilterOptions(gaussian_sigma=0)

    def test_rejects_butterworth_order_below_1(self):
        """Should reject butterworth_order below 1."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            FFTFilterOptions(butterworth_order=0)

    def test_taper_width_bounds(self):
        """Should reject taper_width outside [0, 0.5]."""
        with pytest.raises(ValueError):
            FFTFilterOptions(taper_width=0.6)
        with pytest.raises(ValueError):
            FFTFilterOptions(taper_width=-0.1)

    def test_custom_values(self):
        """Should accept custom values."""
        options = FFTFilterOptions(
            gaussian_sigma=0.5,
            butterworth_order=6,
            spatial_domain=True,
            taper_edges=False,
            taper_width=0.1,
        )

        assert options.gaussian_sigma == pytest.approx(0.5)
        assert options.butterworth_order == 6
        assert options.spatial_domain is True
        assert options.taper_edges is False


class TestCalibrationOptions:
    """Tests for CalibrationOptions configuration."""

    def test_defaults(self):
        """Should have expected default values."""
        options = CalibrationOptions()

        assert options.grid_type == "constant"
        assert options.reference_frame == "IGS20"
        assert options.starting_year == pytest.approx(2014.0)
        assert options.unwrap_error_correction is True
        assert options.window_size_meters == pytest.approx(30000.0)
        assert options.posting_meters == pytest.approx(30.0)
        assert options.downsample_factor == 6
        assert options.downsample_method == "mean"
        assert options.downsample_weighted is False
        assert options.calibration_surface_smoothing_method == "gaussian"
        assert options.calibration_surface_smoothing_sigma is None

    def test_window_size_pixels_computed(self):
        """Should compute window_size_pixels from metres and posting."""
        options = CalibrationOptions(window_size_meters=30000.0, posting_meters=30.0)

        assert options.window_size_pixels == 1000

    def test_custom_grid_type(self):
        """Should accept variable grid type."""
        options = CalibrationOptions(grid_type="variable")

        assert options.grid_type == "variable"

    def test_rejects_invalid_grid_type(self):
        """Should reject unknown grid types."""
        with pytest.raises(ValueError):
            CalibrationOptions(grid_type="unknown")

    def test_rejects_invalid_reference_frame(self):
        """Should reject unknown reference frames."""
        with pytest.raises(ValueError):
            CalibrationOptions(reference_frame="WGS84")

    def test_rejects_zero_downsample_factor(self):
        """Should reject zero downsample factor."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            CalibrationOptions(downsample_factor=0)

    def test_rejects_negative_downsample_factor(self):
        """Should reject negative downsample factor."""
        with pytest.raises(ValueError, match="greater than or equal to 1"):
            CalibrationOptions(downsample_factor=-1)

    def test_rejects_invalid_smoothing_method(self):
        """Should reject unknown smoothing methods."""
        with pytest.raises(ValueError):
            CalibrationOptions(calibration_surface_smoothing_method="cubic")

    def test_all_smoothing_methods(self):
        """Should accept all valid smoothing methods."""
        for method in ("gaussian", "gaussian_fft", "hanning_fft", "savitzky_golay"):
            options = CalibrationOptions(calibration_surface_smoothing_method=method)
            assert options.calibration_surface_smoothing_method == method

    def test_disable_interpolation(self):
        """Should allow setting sigma to zero."""
        options = CalibrationOptions(calibration_surface_smoothing_sigma=0)

        assert options.calibration_surface_smoothing_sigma == 0

    def test_nested_savitzky_golay(self):
        """Should accept nested SavitzkyGolayOptions."""
        options = CalibrationOptions(
            savitzky_golay=SavitzkyGolayOptions(window_length=101, polyorder=5)
        )

        assert options.savitzky_golay.window_length == 101

    def test_nested_fft_filter(self):
        """Should accept nested FFTFilterOptions."""
        options = CalibrationOptions(fft_filter=FFTFilterOptions(gaussian_sigma=0.5))

        assert options.fft_filter.gaussian_sigma == pytest.approx(0.5)

    def test_rejects_non_positive_window_size(self):
        """Should reject non-positive window_size_meters."""
        with pytest.raises(ValueError):
            CalibrationOptions(window_size_meters=0)

    def test_rejects_non_positive_posting(self):
        """Should reject non-positive posting_meters."""
        with pytest.raises(ValueError):
            CalibrationOptions(posting_meters=0)


class TestAlgorithmParameters:
    """Tests for complete AlgorithmParameters."""

    def test_defaults(self):
        """Should create with all default values."""
        params = AlgorithmParameters()

        assert isinstance(params.calibration_options, CalibrationOptions)
        assert params.calibration_options.downsample_factor == 6

    def test_create_default_classmethod(self):
        """Should create defaults via classmethod."""
        params = AlgorithmParameters.create_default()

        assert isinstance(params, AlgorithmParameters)
        assert params.calibration_options.unwrap_error_correction is True

    def test_custom_calibration_options(self):
        """Should accept custom nested configuration."""
        params = AlgorithmParameters(
            calibration_options=CalibrationOptions(
                grid_type="variable",
                downsample_factor=3,
            )
        )

        assert params.calibration_options.grid_type == "variable"
        assert params.calibration_options.downsample_factor == 3

    def test_yaml_round_trip(self, tmp_path: Path):
        """Should serialize and deserialize via YAML."""
        params = AlgorithmParameters.create_default()
        yaml_file = tmp_path / "algorithm_params.yaml"

        params.to_yaml(yaml_file)
        assert yaml_file.exists()

        loaded = AlgorithmParameters.from_yaml(yaml_file)

        assert (
            loaded.calibration_options.calibration_surface_smoothing_method
            == params.calibration_options.calibration_surface_smoothing_method
        )
        assert (
            loaded.calibration_options.downsample_factor
            == params.calibration_options.downsample_factor
        )

    def test_yaml_content_valid(self, tmp_path: Path):
        """Should produce valid YAML output."""
        params = AlgorithmParameters()
        yaml_file = tmp_path / "params.yaml"

        params.to_yaml(yaml_file)

        with open(yaml_file) as f:
            data = yaml.safe_load(f)
        assert "calibration_options" in data

    def test_yaml_with_comments(self, tmp_path: Path):
        """Should generate YAML with field descriptions."""
        params = AlgorithmParameters.create_default()
        yaml_file = tmp_path / "params_commented.yaml"

        params.to_yaml(yaml_file, with_comments=True)

        content = yaml_file.read_text()
        assert "#" in content

    def test_nested_validation_propagates(self):
        """Should propagate validation errors from nested models."""
        with pytest.raises(ValueError):
            AlgorithmParameters(
                calibration_options=CalibrationOptions(downsample_factor=-1)
            )

    def test_modifying_after_creation(self):
        """Should allow modifying fields after creation."""
        params = AlgorithmParameters()

        params.calibration_options.downsample_factor = 20

        assert params.calibration_options.downsample_factor == 20
