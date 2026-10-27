"""Dependency and memory budget of the operational install (Venti plan T26.4, T26.5)."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def budget_mod():
    spec = importlib.util.spec_from_file_location(
        "check_env_budget", ROOT / "scripts" / "check_env_budget.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_check_classifies_violations_tolerated_and_stale(budget_mod):
    budget = {
        "forbidden": ["dask", "zarr", "matplotlib", "numba"],
        "known_exceptions": {"matplotlib": "browse", "numba": "old reason"},
    }
    violations, tolerated, stale = budget_mod.check(
        {"dask", "matplotlib", "numpy"}, budget
    )
    assert violations == ["dask"]
    assert tolerated == ["matplotlib"]
    assert stale == ["numba"]


def test_core_install_has_no_dask_or_zarr(budget_mod):
    """Regression: opera-utils[disp] in the core brought dask, zarr and s3fs for
    one reference-pixel function; the rule is ported and the extra moved to
    `download`."""
    closure, _missing = budget_mod.installed_closure("cal-disp")
    assert "numpy" in closure and "venti" in closure
    assert "dask" not in closure
    assert "zarr" not in closure
    assert "asf_search" not in closure  # download extra only


def test_cli_passes_on_the_core_install(budget_mod, capsys):
    assert budget_mod.main(["--closure", "cal-disp"]) == 0
    out = capsys.readouterr().out
    assert "tolerated  matplotlib" in out and "OK" in out


def test_report_mode_reads_a_pip_resolution(budget_mod, tmp_path, capsys):
    report = {
        "install": [{"metadata": {"name": n}} for n in ("numpy", "Dask", "matplotlib")]
    }
    path = tmp_path / "r.json"
    path.write_text(json.dumps(report))
    assert budget_mod.main(["--report", str(path)]) == 1
    assert "FORBIDDEN  dask" in capsys.readouterr().out


def test_peak_memory_against_the_budget(budget_mod, tmp_path, capsys):
    log = tmp_path / "run.log"
    log.write_text(
        "...\nMaximum memory usage: 9.10 GB\n...\nMaximum memory usage: 3.95 GB\n"
    )
    assert budget_mod.peak_memory_gb(log.read_text()) == 3.95
    assert budget_mod.main(["--peak-memory-log", str(log)]) == 0
    log.write_text("Maximum memory usage: 6.10 GB\n")  # above the 6 GB budget
    assert budget_mod.main(["--peak-memory-log", str(log)]) == 1
    assert "FAILED" in capsys.readouterr().out
    with pytest.raises(ValueError, match="Maximum memory usage"):
        budget_mod.peak_memory_gb("nothing here")
