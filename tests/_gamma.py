"""The gamma (0.3) parameters for tests of the gamma route.

Since v0.5 (plan T49) the defaults are the v0.5 algorithm; a test of the gamma
route says so with these parameters instead of relying on the defaults.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from cal_disp.config._algorithm import AlgorithmParameters

GAMMA_FILE = (
    Path(__file__).resolve().parents[1] / "configs" / "algorithm_parameters_gamma.yaml"
)


def gamma_parameters(overrides: dict[str, Any] | None = None) -> AlgorithmParameters:
    """`configs/algorithm_parameters_gamma.yaml`, with optional overrides."""
    return AlgorithmParameters.from_yaml(GAMMA_FILE).with_overrides(overrides or {})
