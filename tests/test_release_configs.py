"""Release parameter files (Venti plan T49): gamma reproducible, v0.5 = T46 set."""

from pathlib import Path

import yaml

from cal_disp.config._algorithm import AlgorithmParameters

ROOT = Path(__file__).resolve().parents[1]


def _golden_03_config() -> dict:
    """The 0.3 golden's parameters (frozen in tests/ when T49 moved to v0.5)."""
    data = yaml.safe_load(
        (ROOT / "tests" / "golden_0.3_algorithm_parameters.yaml").read_text()
    )
    # the 0.3 golden relied on the then-defaults of the schema-v2 groups: Venti's
    from venti.workflow.config import (
        GnssOptions,
        SurfaceOptions,
        UncertaintyOptions,
        WeightOptions,
    )

    for key, model in (
        ("surface", SurfaceOptions),
        ("weights", WeightOptions),
        ("gnss", GnssOptions),
        ("uncertainty", UncertaintyOptions),
    ):
        data["calibration_options"][key] = model().model_dump()
    return data


def _load(name: str) -> AlgorithmParameters:
    return AlgorithmParameters.model_validate(
        yaml.safe_load((ROOT / "configs" / name).read_text())
    )


def test_gamma_file_is_the_golden_configuration():
    """Decision 0001: the gamma file reproduces the 0.3 golden's parameters."""
    golden = AlgorithmParameters.model_validate(_golden_03_config())
    assert _load("algorithm_parameters_gamma.yaml") == golden
    cal = golden.calibration_options
    assert cal.window_size_meters == 600_000.0
    assert cal.apply_tropo_correction is False
    assert cal.unwrap_error_correction is False
    assert cal.surface.method == "windowed_plane"


def test_v05_file_is_the_t46_benchmark_set():
    cal = _load("algorithm_parameters_v05.yaml").calibration_options
    assert cal.surface.method == "loclin"
    assert cal.surface.cutoff_wavelength_meters == 50_000.0
    assert cal.surface.two_pass and cal.surface.fill_gaps
    assert (cal.weights.coherence_power, cal.weights.robust) == (8.0, True)
    assert cal.gnss.buffer_meters == 50_000.0
    assert cal.gnss.exclude_defo_nodes and cal.gnss.reinterpolate_excluded
    # per-frame settings come from the frame table, not this file
    assert cal.surface.regions == "none"
    assert cal.surface.fit_mask is None
    # owner 2026-10-14: sigma_CAL wherever the calibration has a value
    assert cal.uncertainty.extrapolate_sigma is True
    # owner 2026-10-09: off until a trade study confirms it
    assert cal.unwrap_error_correction is False


def test_model_defaults_are_the_v05_set():
    """T49.1: the defaults are the v0.5 algorithm; gamma needs its own file."""
    assert AlgorithmParameters() == _load("algorithm_parameters_v05.yaml")
    assert AlgorithmParameters() != _load("algorithm_parameters_gamma.yaml")


def test_golden_build_writes_the_v05_set_and_the_frame_table():
    """T49.1: the golden is built with the v0.5 file and production's ancillaries."""
    script = (ROOT / "scripts" / "build_golden_output.sh").read_text()
    assert "configs/algorithm_parameters_v05.yaml" in script
    for flag in ("--algorithm-overrides", "--defo-area-db", "--event-db"):
        assert flag in script


def test_shipped_config_is_the_defaults():
    """configs/algorithm_parameters.yaml documents the defaults (the v0.5 set)."""
    assert _load("algorithm_parameters.yaml") == AlgorithmParameters()
