# Placeholder for the algorithm options
# _algorithm.py
from __future__ import annotations

from enum import Enum
from typing import Optional

from pydantic import Field

from ._yaml import YamlModel


class CalibrationMethod(str, Enum):
    """Low pass filter methods for calibration plane fitting."""

    GAUSSIAN_FFT = "gaussian_fft"
    HANNING_FFT = "hanning_fft"
    SAVITSKY_GOLEY = "savitsky_goley"


class UnwrapCorrection(YamlModel):
    """Unwrapping error correction configuration.

    Attributes
    ----------
    run_unwrap_correction : bool
        Whether to run the unwrapping error correction before calibration.

    """

    run_unwrap_correction: bool = Field(
        default=True,
        description=(
            "Whether to run the unwrapping error correction before calibration."
        ),
    )


class SavitskyGoleyOptions(YamlModel):
    """Savitsky-Golay filter options for calibration.

    Attributes
    ----------
    window_x_size : int | None
        Moving window size for plane fitting along x-axis (columns), in pixels.
    window_y_size : int | None
        Moving window size for plane fitting along y-axis (rows), in pixels.
    window_overlap_x_size : int | None
        Overlap for x-axis windows, in pixels.
    window_overlap_y_size : int | None
        Overlap for y-axis windows, in pixels.
    window_extend_x_size : int | None
        Moving x-axis window buffer to extend kernel, in pixels.
    window_extend_y_size : int | None
        Moving y-axis window buffer to extend kernel, in pixels.

    """

    window_x_size: Optional[int] = Field(
        default=None,
        description=(
            "Moving window size for plane fitting along x-axis (columns), in pixels."
        ),
        ge=1,
    )

    window_y_size: Optional[int] = Field(
        default=None,
        description=(
            "Moving window size for plane fitting along y-axis (rows), in pixels."
        ),
        ge=1,
    )

    window_overlap_x_size: Optional[int] = Field(
        default=None,
        description="Overlap for x-axis windows, in pixels.",
        ge=0,
    )

    window_overlap_y_size: Optional[int] = Field(
        default=None,
        description="Overlap for y-axis windows, in pixels.",
        ge=0,
    )

    window_extend_x_size: Optional[int] = Field(
        default=None,
        description="Moving x-axis window buffer to extend kernel, in pixels.",
        ge=0,
    )

    window_extend_y_size: Optional[int] = Field(
        default=None,
        description="Moving y-axis window buffer to extend kernel, in pixels.",
        ge=0,
    )


class CalibrationOptions(YamlModel):
    """Calibration algorithm configuration options.

    Attributes
    ----------
    cal_method : CalibrationMethod
        Low pass method to get a long wavelength calibration plane.
    run_interpolation : bool
        Whether to run gap filling interpolation step before lowpass.
    run_downsample : bool
        Whether to downsample DISP before calibration.
    downsample_factor : int
        Downsampling factor for DISP.

    """

    cal_method: CalibrationMethod = Field(
        default=CalibrationMethod.SAVITSKY_GOLEY,
        description=(
            "Low pass method to get a long wavelength calibration plane. "
            "Options: gaussian_fft, hanning_fft, savitsky_goley."
        ),
    )

    run_interpolation: bool = Field(
        default=True,
        description="Whether to run gap filling interpolation step before lowpass.",
    )

    run_downsample: bool = Field(
        default=True,
        description="Whether to downsample DISP before calibration.",
    )

    downsample_factor: int = Field(
        default=10,
        description="Downsampling factor for DISP.",
        ge=1,
    )


class AlgorithmParameters(YamlModel):
    """CAL-DISP Algorithm Parameters.

    Configuration for the calibration algorithm.

    Attributes
    ----------
    unwrap_correction : UnwrapCorrection
        Unwrapping error correction configuration.
    calibration_options : CalibrationOptions
        Calibration algorithm configuration options.
    savitsky_goley_options : SavitskyGoleyOptions
        Savitsky-Golay filter specific options.

    Examples
    --------
    >>> # Create with defaults
    >>> params = AlgorithmParameters()
    >>>
    >>> # Load from YAML
    >>> params = AlgorithmParameters.from_yaml("algorithm_params.yaml")
    >>>
    >>> # Modify and save
    >>> params.calibration_options.downsample_factor = 5
    >>> params.to_yaml("updated_params.yaml")

    """

    unwrap_correction: UnwrapCorrection = Field(
        default_factory=UnwrapCorrection,
        description="Unwrapping error correction configuration.",
    )

    calibration_options: CalibrationOptions = Field(
        default_factory=CalibrationOptions,
        description="Calibration algorithm configuration options.",
    )

    savitsky_goley_options: SavitskyGoleyOptions = Field(
        default_factory=SavitskyGoleyOptions,
        description="Savitsky-Golay filter specific options.",
    )

    @classmethod
    def create_default(cls) -> "AlgorithmParameters":
        """Create algorithm parameters with default values.

        Returns
        -------
        AlgorithmParameters
            Default configuration.

        """
        return cls()

    @classmethod
    def create_example(cls) -> "AlgorithmParameters":
        """Create example algorithm parameters with typical values.

        Returns
        -------
        AlgorithmParameters
            Example configuration with all options specified.

        """
        return cls(
            unwrap_correction=UnwrapCorrection(run_unwrap_correction=True),
            calibration_options=CalibrationOptions(
                cal_method=CalibrationMethod.SAVITSKY_GOLEY,
                run_interpolation=True,
                run_downsample=True,
                downsample_factor=10,
            ),
            savitsky_goley_options=SavitskyGoleyOptions(
                window_x_size=100,
                window_y_size=100,
                window_overlap_x_size=50,
                window_overlap_y_size=50,
                window_extend_x_size=10,
                window_extend_y_size=10,
            ),
        )
