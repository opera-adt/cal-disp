"""Algorithm parameter configuration for CAL-DISP."""

from __future__ import annotations

from typing import Literal, Optional

from pydantic import Field, computed_field

from ._yaml import YamlModel


class SavitzkyGolayOptions(YamlModel):
    """Savitzky-Golay filter options for calibration surface smoothing.

    Attributes
    ----------
    window_length : int
        Filter window length in pixels (must be odd, >= 3).
    polyorder : int
        Polynomial order for fitting (must be less than window_length).
    deriv : int
        Derivative order (0 = smoothing only).

    """

    window_length: int = Field(
        default=51,
        ge=3,
        description="Filter window length in pixels (must be odd).",
    )
    polyorder: int = Field(
        default=3,
        ge=0,
        description="Polynomial order for fitting.",
    )
    deriv: int = Field(
        default=0,
        ge=0,
        description="Derivative order (0 = smoothing only).",
    )


class FFTFilterOptions(YamlModel):
    """FFT-based filter options for calibration surface smoothing.

    Attributes
    ----------
    gaussian_sigma : float
        Gaussian filter standard deviation in wavelength units.
    butterworth_order : int
        Butterworth filter order.
    spatial_domain : bool
        Apply filter in spatial domain instead of frequency domain.
    taper_edges : bool
        Taper image edges to reduce ringing artefacts.
    taper_width : float
        Taper width as a fraction of the image dimensions [0, 0.5].

    """

    gaussian_sigma: float = Field(
        default=0.1,
        gt=0,
        description="Gaussian filter standard deviation in wavelength units.",
    )
    butterworth_order: int = Field(
        default=4,
        ge=1,
        description="Butterworth filter order.",
    )
    spatial_domain: bool = Field(
        default=False,
        description="Apply filter in spatial domain instead of frequency domain.",
    )
    taper_edges: bool = Field(
        default=True,
        description="Taper image edges to reduce ringing artefacts.",
    )
    taper_width: float = Field(
        default=0.05,
        ge=0.0,
        le=0.5,
        description="Taper width as a fraction of the image dimensions.",
    )


class CalibrationOptions(YamlModel):
    """Calibration algorithm configuration.

    Controls GNSS grid type, windowed plane-fitting window size, downsampling,
    unwrap-error correction, and post-assembly smoothing of the calibration
    surface.  These parameters map directly to the Venti ``GNSSReference`` and
    ``SpatialProcessor`` API.

    Attributes
    ----------
    grid_type : {'constant', 'variable'}
        GNSS model type.  ``'constant'`` interpolates a single velocity field
        scaled to the acquisition interval; ``'variable'`` fetches epoch-specific
        GNSS displacements for each (ref_date, sec_date) pair.
    reference_frame : str
        GNSS reference frame, ``'IGS20'`` or ``'IGS14'``.
    starting_year : float
        Earliest observation year used to estimate station velocities
        (``constant`` grid type only).
    unwrap_error_correction : bool
        Apply watershed-based unwrap-error correction before fitting.
    window_size_meters : float
        Side length of the moving-window used for polynomial plane fitting,
        in metres.
    posting_meters : float
        Input DISP pixel spacing in metres (30 m for DISP-S1).
    downsample_factor : int
        Integer downsampling factor applied before surface fitting.  Set to 1
        to disable.  The calibration surface is upsampled back to full
        resolution after fitting.
    downsample_method : {'mean', 'median'}
        Aggregation method used when downsampling.
    calibration_surface_smoothing_method : {'gaussian', 'gaussian_fft', 'hanning_fft', 'savitzky_golay'}
        Post-assembly low-pass filter applied to the stitched calibration
        surface to suppress window-boundary artefacts.
    calibration_surface_smoothing_sigma : float or None
        Sigma (pixels) for ``'gaussian'`` and FFT smoothing methods.
        ``None`` (default) auto-selects ``window_size_pixels / 8``; ``0``
        disables smoothing.  Ignored when method is ``'savitzky_golay'``.
    savitzky_golay : SavitzkyGolayOptions
        Savitzky-Golay filter parameters, used when
        ``calibration_surface_smoothing_method = 'savitzky_golay'``.
    fft_filter : FFTFilterOptions
        FFT filter parameters, used when
        ``calibration_surface_smoothing_method`` is an FFT variant.

    """

    grid_type: Literal["constant", "variable"] = Field(
        default="constant",
        description=(
            "GNSS model type: 'constant' uses a velocity field scaled to the "
            "acquisition interval; 'variable' fetches epoch-specific displacements."
        ),
    )

    reference_frame: Literal["IGS20", "IGS14"] = Field(
        default="IGS20",
        description="GNSS reference frame used for station timeseries.",
    )

    starting_year: float = Field(
        default=2014.0,
        description=(
            "Earliest observation year included in velocity estimation "
            "(constant grid type only)."
        ),
    )

    unwrap_error_correction: bool = Field(
        default=True,
        description=(
            "Apply watershed-based unwrap-error correction to the displacement "
            "field before fitting the calibration surface."
        ),
    )

    window_size_meters: float = Field(
        default=30000.0,
        gt=0,
        description=(
            "Side length of the moving window used for polynomial plane fitting, "
            "in metres."
        ),
    )

    posting_meters: float = Field(
        default=30.0,
        gt=0,
        description="Input DISP pixel spacing in metres (30 m for DISP-S1).",
    )

    downsample_factor: int = Field(
        default=6,
        ge=1,
        description=(
            "Integer downsampling factor applied before surface fitting. "
            "Set to 1 to disable.  The calibration surface is upsampled back "
            "to full resolution after fitting."
        ),
    )

    downsample_method: Literal["mean", "median"] = Field(
        default="mean",
        description="Pixel aggregation method used during downsampling.",
    )

    downsample_weighted: bool = Field(
        default=False,
        description=(
            "Weight downsampling by temporal coherence. "
            "When True, the temporal_coherence layer from the DISP product is "
            "used as per-pixel weights during aggregation.  Falls back to "
            "unweighted downsampling if the layer is unavailable."
        ),
    )

    calibration_surface_smoothing_method: Literal[
        "gaussian", "gaussian_fft", "hanning_fft", "savitzky_golay"
    ] = Field(
        default="gaussian",
        description=(
            "Post-assembly low-pass filter applied to the stitched calibration "
            "surface.  Options: 'gaussian' (spatial-domain, default), "
            "'gaussian_fft', 'hanning_fft', or 'savitzky_golay'."
        ),
    )

    calibration_surface_smoothing_sigma: Optional[float] = Field(
        default=None,
        ge=0,
        description=(
            "Sigma (pixels) for gaussian and FFT smoothing methods.  "
            "None auto-selects window_size_pixels/8; 0 disables smoothing.  "
            "Ignored when calibration_surface_smoothing_method='savitzky_golay'."
        ),
    )

    savitzky_golay: SavitzkyGolayOptions = Field(
        default_factory=SavitzkyGolayOptions,
        description=(
            "Savitzky-Golay filter parameters.  Active when "
            "calibration_surface_smoothing_method='savitzky_golay'."
        ),
    )

    fft_filter: FFTFilterOptions = Field(
        default_factory=FFTFilterOptions,
        description=(
            "FFT filter parameters.  Active when "
            "calibration_surface_smoothing_method is an FFT variant."
        ),
    )

    @computed_field  # type: ignore[misc]
    @property
    def window_size_pixels(self) -> int:
        """Window size in pixels derived from metres and posting."""
        return max(1, int(self.window_size_meters / self.posting_meters))


class AlgorithmParameters(YamlModel):
    """CAL-DISP Algorithm Parameters.

    Top-level container for all algorithm configuration.  Load from YAML with
    ``AlgorithmParameters.from_yaml(path)``; serialise with ``.to_yaml(path)``.

    Attributes
    ----------
    calibration_options : CalibrationOptions
        Controls GNSS grid type, windowed surface fitting, and downsampling.

    Examples
    --------
    >>> params = AlgorithmParameters()
    >>> params.calibration_options.window_size_meters
    30000.0
    >>> params = AlgorithmParameters.from_yaml("algorithm_parameters.yaml")

    """

    calibration_options: CalibrationOptions = Field(
        default_factory=CalibrationOptions,
        description="GNSS calibration and windowed surface-fitting options.",
    )

    @classmethod
    def create_default(cls) -> "AlgorithmParameters":
        """Create algorithm parameters with default values."""
        return cls()
