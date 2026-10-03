"""Algorithm parameter configuration for CAL-DISP."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Literal, Optional

from pydantic import Field, computed_field

from ._yaml import YamlModel, _computed_fields_exclude


class SavitzkyGolayOptions(YamlModel):
    """Savitzky-Golay filter options for calibration surface smoothing.

    Attributes
    ----------
    window_length : int
        Filter window length in pixels (must be odd, >= 3).
    polyorder : int
        Polynomial order for fitting (must be less than window_length).

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


class CalibrationOptions(YamlModel):
    """Calibration algorithm configuration.

    Controls GNSS grid type, windowed plane-fitting window size, downsampling,
    unwrap-error correction, and post-assembly smoothing of the calibration
    surface.  `to_venti` converts them to Venti's ``CalibrationOptions`` for
    ``venti.estimate_calibration_surface``; the downsampling options are
    passed to it as arguments.

    Attributes
    ----------
    grid_type : {'constant', 'variable'}
        UNR grid product.  ``'constant'``: UNR's precomputed linear rates times
        the acquisition interval; ``'variable'``: the difference of the daily
        positions nearest the reference and secondary dates.
    reference_frame : str
        GNSS reference frame, ``'IGS20'`` or ``'IGS14'``.
    unwrap_error_correction : bool
        Apply Venti's region-offset unwrap-error correction before fitting.
        Off by default: Venti segments the *valid mask* into connected
        components (not phase discontinuities) and shifts every component
        whose median differs from the first one by more than half a cycle
        (λ/2 = 27.7 mm for Sentinel-1), which quantises genuine
        long-wavelength signal into λ/2 steps before the plane fit.
    apply_tropo_correction : bool
        Apply the tropospheric correction when the runconfig lists tropo
        files; ``False`` ignores them.
    apply_solid_earth_tide_correction : bool
        Remove the product's ``/corrections/solid_earth_tide`` before the fit
        (UNR GNSS already has it removed) and add it back to the surface.
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
    downsample_weighted : bool
        Weight downsampling by the product's ``temporal_coherence``.
    event_mask_buffer_pixels : int
        Pixels by which event / deformation masks are dilated.
    residual_outlier_mad_threshold : float or None
        Without an event mask, fill pixels whose InSAR - GNSS residual exceeds
        this many robust sigmas before the fit.  ``None`` disables.
    residual_region_mad_threshold : float or None
        Without an event mask, fill coherent residual regions above this many
        robust sigmas.  ``None`` disables.
    residual_region_min_pixels : int
        Minimum size of such a region.
    mask_fit_residual_outliers : bool
        Trim the most extreme 15% / 85% residual quantiles inside each window.
    weight_fit_by_gnss_uncertainty : bool
        Weight the windowed fit by the inverse GNSS LOS variance.
    calibration_surface_smoothing_method : {'gaussian', 'gaussian_fft',
        'hanning_fft', 'savitzky_golay'}
        Post-assembly low-pass filter applied to the stitched calibration
        surface to suppress window-boundary artefacts.
    calibration_surface_smoothing_sigma : float or None
        Sigma (pixels) for ``'gaussian'`` and FFT smoothing methods.
        ``None`` (default) auto-selects ``window_size_pixels / 8``; ``0``
        disables smoothing.  Ignored when method is ``'savitzky_golay'``.
    savitzky_golay : SavitzkyGolayOptions
        Savitzky-Golay filter parameters, used when
        ``calibration_surface_smoothing_method = 'savitzky_golay'``.

    """

    grid_type: Literal["constant", "variable"] = Field(
        default="constant",
        description=(
            "UNR grid product: 'constant' uses UNR's precomputed rates scaled to "
            "the acquisition interval; 'variable' uses the positions nearest the "
            "two dates. Must match the staged UNR data (unr_grid_type)."
        ),
    )

    reference_frame: Literal["IGS20", "IGS14"] = Field(
        default="IGS20",
        description="GNSS reference frame used for station timeseries.",
    )

    unwrap_error_correction: bool = Field(
        default=False,
        description=(
            "Apply Venti's region-offset unwrap-error correction to the "
            "displacement field before fitting the calibration surface. "
            "Disabled by default: Venti segments the valid-data mask into "
            "connected components (islands), not into regions bounded by phase "
            "discontinuities, and shifts every island whose median differs from "
            "the first island by more than half an unwrapping cycle (lambda/2 = "
            "27.7 mm for Sentinel-1). On a real frame this quantises the "
            "long-wavelength signal the calibration is meant to fit into "
            "lambda/2 steps. Enable only once the segmentation is driven by "
            "phase discontinuities (or connected components) and validated "
            "against products with known unwrapping errors."
        ),
    )

    apply_tropo_correction: bool = Field(
        default=True,
        description=(
            "Apply the tropospheric correction when the runconfig lists tropo "
            "files (removed before the fit and added back to the surface). "
            "false ignores the files."
        ),
    )

    apply_solid_earth_tide_correction: bool = Field(
        default=True,
        description=(
            "Remove the product's /corrections/solid_earth_tide before the fit "
            "(UNR GNSS already has it removed) and add it back to the surface. "
            "A product without the layer is calibrated without it (warning)."
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
            "used as per-pixel weights during aggregation; an error is raised "
            "if the layer is missing."
        ),
    )

    event_mask_buffer_pixels: int = Field(
        default=0,
        ge=0,
        description="Pixels by which event / deformation masks are dilated.",
    )

    residual_outlier_mad_threshold: Optional[float] = Field(
        default=None,
        gt=0,
        description=(
            "Only without an event mask: fill pixels whose InSAR - GNSS residual "
            "exceeds this many robust (1.4826 * MAD) sigmas before the fit. "
            "null disables."
        ),
    )

    residual_region_mad_threshold: Optional[float] = Field(
        default=None,
        gt=0,
        description=(
            "Only without an event mask: fill coherent residual regions above "
            "this many robust sigmas before the fit. null disables."
        ),
    )

    residual_region_min_pixels: int = Field(
        default=20,
        ge=1,
        description="Minimum size (downsampled pixels) of a residual region.",
    )

    mask_fit_residual_outliers: bool = Field(
        default=True,
        description=(
            "Trim the most extreme 15% / 85% InSAR - GNSS residual quantiles "
            "inside each fit window."
        ),
    )

    weight_fit_by_gnss_uncertainty: bool = Field(
        default=False,
        description="Weight the windowed fit by the inverse GNSS LOS variance.",
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

    @computed_field  # type: ignore[misc]
    @property
    def window_size_pixels(self) -> int:
        """Window size in pixels derived from metres and posting."""
        return max(1, int(self.window_size_meters / self.posting_meters))

    def to_venti(self):
        """Venti ``CalibrationOptions`` with these settings.

        Options Venti does not take from its config (downsampling) are
        passed to ``estimate_calibration_surface`` separately.
        """
        from venti.workflow.config import CalibrationOptions as VentiOptions

        return VentiOptions(
            **self.model_dump(include=VENTI_OPTIONS - {"savitzky_golay"}),
            savitzky_golay=self.savitzky_golay.model_dump(),
        )


# CalibrationOptions fields passed to Venti's CalibrationOptions
VENTI_OPTIONS = {
    "grid_type",
    "reference_frame",
    "unwrap_error_correction",
    "apply_tropo_correction",
    "apply_solid_earth_tide_correction",
    "window_size_meters",
    "posting_meters",
    "event_mask_buffer_pixels",
    "residual_outlier_mad_threshold",
    "residual_region_mad_threshold",
    "residual_region_min_pixels",
    "mask_fit_residual_outliers",
    "weight_fit_by_gnss_uncertainty",
    "calibration_surface_smoothing_method",
    "calibration_surface_smoothing_sigma",
    "savitzky_golay",
}
# Venti options cal-disp does not expose: each run has its own GNSS cache.
VENTI_OPTIONS_NOT_EXPOSED = {"recompute_gnss"}


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

    def with_overrides(self, overrides: Mapping[str, Any]) -> "AlgorithmParameters":
        """Return a validated copy with `overrides` applied.

        Parameters
        ----------
        overrides : Mapping[str, Any]
            Frame-specific values, e.g. one entry of the
            ``algorithm_parameters_overrides_json`` file. Keys are either
            option groups (``{"calibration_options": {"downsample_factor": 3}}``)
            or, as a shorthand, options of ``calibration_options``
            (``{"downsample_factor": 3}``). Nested options are merged.

        Returns
        -------
        AlgorithmParameters
            New instance; ``self`` is returned unchanged for empty overrides.

        Raises
        ------
        pydantic.ValidationError
            For an unknown option or an invalid value.

        """
        if not overrides:
            return self

        def _merge(base: dict[str, Any], new: Mapping[str, Any]) -> dict[str, Any]:
            for key, value in new.items():
                if isinstance(value, Mapping) and isinstance(base.get(key), dict):
                    _merge(base[key], value)
                else:
                    base[key] = value
            return base

        # Without computed fields: they would be rejected as extra inputs
        data = self.model_dump(exclude=_computed_fields_exclude(self))
        groups = set(type(self).model_fields)
        for key, value in overrides.items():
            if key in groups:
                _merge(data, {key: value})
            else:
                _merge(data["calibration_options"], {key: value})
        return type(self).model_validate(data)
