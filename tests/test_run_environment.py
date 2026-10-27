"""``cal-disp run`` as a PGE would: separate process, restricted directories."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from click.testing import CliRunner

from cal_disp.cli import cli_app


def _write_runconfig(tmp_path, name, disp, los, dem, unr, algo) -> Path:
    """Write ``<name>/work/runconfig.yaml`` (the CLI writes it to the work dir)."""
    lookup_file, tenv8_dir = unr
    run_dir = tmp_path / name
    work_dir, output_dir = run_dir / "work", run_dir / "output"
    work_dir.mkdir(parents=True)
    output_dir.mkdir()
    config_file = work_dir / "runconfig.yaml"
    result = CliRunner().invoke(
        cli_app,
        # fmt: off
        [
            "config", "-d", str(disp), "-ul", str(lookup_file),
            "-ud", str(tenv8_dir), "-uv", "0.3", "-ut", "constant",
            "--los-file", str(los), "--dem-file", str(dem), "-a", str(algo),
            "-c", str(config_file), "--frame-id", "8882",
            "-o", str(output_dir), "--work-dir", str(work_dir),
        ],
        # fmt: on
    )
    assert result.exit_code == 0, result.output
    return config_file


def _run(config_file: Path, cwd: Path, env: dict[str, str]) -> Path:
    proc = subprocess.run(
        [sys.executable, "-c", "from cal_disp.cli import cli_app; cli_app()"]
        + ["run", str(config_file)],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    (output,) = (config_file.parent.parent / "output").glob("*.nc")
    return output


@pytest.mark.slow  # two full `cal-disp run` subprocesses
@pytest.mark.skipif(
    hasattr(os, "geteuid") and os.geteuid() == 0, reason="root ignores permissions"
)
def test_run_without_writable_cwd_or_system_temp(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
):
    """GDAL/tempfile work files must go to the scratch dir, not cwd or /tmp.

    Reproduces a container where the working directory (e.g. /data/work) and
    the system temp dirs are not writable: GDAL then failed silently inside
    Venti's gap filling. The output must match an unrestricted run exactly.
    """
    inputs = (
        sample_disp_product,
        sample_static_los,
        sample_static_dem,
        sample_unr_data,
        sample_algorithm_params,
    )
    env = {k: v for k, v in os.environ.items() if k != "CPL_TMPDIR"}
    env["OMP_NUM_THREADS"] = "1"
    reference = _run(_write_runconfig(tmp_path, "open", *inputs), tmp_path, env)

    locked = tmp_path / "locked"
    locked.mkdir()
    locked.chmod(0o555)
    restricted_env = {**env, "TMPDIR": str(locked), "TEMP": str(locked)}
    restricted_env["TMP"] = str(locked)
    config_file = _write_runconfig(tmp_path, "restricted", *inputs)
    try:
        output = _run(config_file, locked, restricted_env)
    finally:
        locked.chmod(0o755)

    assert list(locked.iterdir()) == []
    with xr.open_dataset(reference) as ref, xr.open_dataset(output) as out:
        np.testing.assert_array_equal(out["calibration"], ref["calibration"])
        assert np.isfinite(out["calibration"]).any()
    # The run's private temp folder under <scratch>/tmp is removed afterwards
    assert not list(config_file.parent.glob("tmp/run_*"))

    # Venti's messages reach the run's log file too
    log_text = (config_file.parent.parent / "output" / "cal_disp.log").read_text()
    assert " - cal_disp." in log_text
    assert " - venti." in log_text
