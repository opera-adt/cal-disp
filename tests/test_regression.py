"""Regression tests for cal-disp using conftest fixtures.

Tests the complete workflow: config generation -> run -> validate outputs.
Uses existing fixtures from conftest.py to create sample inputs.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr
from click.testing import CliRunner

# Import CLI app (it's named cli_app in the module)
try:
    from cal_disp.cli import cli_app as cli
except ImportError:
    cli = None


@pytest.mark.integration
def test_full_workflow(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
    cli_runner: CliRunner,
) -> None:
    """Test complete cal-disp workflow from config to output.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory.
    sample_disp_product : Path
        Mock DISP-S1 NetCDF file.
    sample_static_los : Path
        Mock LOS GeoTIFF with 3 bands.
    sample_static_dem : Path
        Mock DEM GeoTIFF.
    sample_unr_data : tuple[Path, Path]
        (lookup_file, tenv8_dir) for UNR GPS data.
    sample_algorithm_params : Path
        Sample algorithm parameters YAML file.
    cli_runner : CliRunner
        Click CLI test runner.

    """
    if cli is None:
        pytest.skip("Could not import CLI module")

    lookup_file, tenv8_dir = sample_unr_data
    output_dir = tmp_path / "output"
    work_dir = tmp_path / "work"
    output_dir.mkdir()
    work_dir.mkdir()

    # Generate config
    config_file = work_dir / "runconfig.yaml"
    result = cli_runner.invoke(
        cli,
        [
            "config",
            "-d",
            str(sample_disp_product),
            "-ul",
            str(lookup_file),
            "-ud",
            str(tenv8_dir),
            "-uv",
            "0.2",
            "-ut",
            "variable",
            "--los-file",
            str(sample_static_los),
            "--dem-file",
            str(sample_static_dem),
            "-a",
            str(sample_algorithm_params),
            "-c",
            str(config_file),
            "--frame-id",
            "8882",
            "-o",
            str(output_dir),
            "--work-dir",
            str(work_dir),
        ],
    )

    assert result.exit_code == 0, f"Config failed: {result.output}"
    assert config_file.exists(), "Config file not created"

    # Run calibration
    result = cli_runner.invoke(cli, ["run", str(config_file)])
    assert result.exit_code == 0, f"Run failed: {result.output}"

    # Validate outputs exist
    output_files = list(output_dir.glob("*.nc"))
    assert len(output_files) > 0, "No output NetCDF files created"

    # Validate output structure
    for output_file in output_files:
        _validate_output(output_file)


@pytest.mark.integration
def test_workflow_with_corrections(
    tmp_path: Path,
    sample_disp_product_with_corrections: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
    cli_runner: CliRunner,
) -> None:
    """Test workflow with DISP product that has corrections group.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory.
    sample_disp_product_with_corrections : Path
        Mock DISP-S1 file with corrections group.
    sample_static_los : Path
        Mock LOS GeoTIFF.
    sample_static_dem : Path
        Mock DEM GeoTIFF.
    sample_unr_data : tuple[Path, Path]
        (lookup_file, tenv8_dir).
    sample_algorithm_params : Path
        Sample algorithm parameters YAML file.
    cli_runner : CliRunner
        Click CLI test runner.

    """
    if cli is None:
        pytest.skip("Could not import CLI module")

    lookup_file, tenv8_dir = sample_unr_data
    output_dir = tmp_path / "output"
    work_dir = tmp_path / "work"
    output_dir.mkdir()
    work_dir.mkdir()

    config_file = work_dir / "runconfig.yaml"
    result = cli_runner.invoke(
        cli,
        [
            "config",
            "-d",
            str(sample_disp_product_with_corrections),
            "-ul",
            str(lookup_file),
            "-ud",
            str(tenv8_dir),
            "-uv",
            "0.2",
            "-ut",
            "variable",
            "--los-file",
            str(sample_static_los),
            "--dem-file",
            str(sample_static_dem),
            "-a",
            str(sample_algorithm_params),
            "-c",
            str(config_file),
            "--frame-id",
            "8882",
            "-o",
            str(output_dir),
            "--work-dir",
            str(work_dir),
        ],
    )

    assert result.exit_code == 0
    assert config_file.exists()

    result = cli_runner.invoke(cli, ["run", str(config_file)])
    assert result.exit_code == 0

    output_files = list(output_dir.glob("*.nc"))
    assert len(output_files) > 0


@pytest.mark.integration
def test_workflow_idempotency(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
    cli_runner: CliRunner,
) -> None:
    """Test that running workflow twice produces identical results.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory.
    sample_disp_product : Path
        Mock DISP file.
    sample_static_los : Path
        Mock LOS file.
    sample_static_dem : Path
        Mock DEM file.
    sample_unr_data : tuple[Path, Path]
        (lookup_file, tenv8_dir).
    sample_algorithm_params : Path
        Sample algorithm parameters YAML file.
    cli_runner : CliRunner
        Click CLI test runner.

    """
    if cli is None:
        pytest.skip("Could not import CLI module")

    lookup_file, tenv8_dir = sample_unr_data

    # Run workflow twice
    outputs = []
    for run_id in [1, 2]:
        output_dir = tmp_path / f"run{run_id}" / "output"
        work_dir = tmp_path / f"run{run_id}" / "work"
        output_dir.mkdir(parents=True)
        work_dir.mkdir(parents=True)

        config_file = work_dir / "runconfig.yaml"
        result = cli_runner.invoke(
            cli,
            [
                "config",
                "-d",
                str(sample_disp_product),
                "-ul",
                str(lookup_file),
                "-ud",
                str(tenv8_dir),
                "-uv",
                "0.2",
                "-ut",
                "variable",
                "--los-file",
                str(sample_static_los),
                "--dem-file",
                str(sample_static_dem),
                "-a",
                str(sample_algorithm_params),
                "-c",
                str(config_file),
                "--frame-id",
                "8882",
                "-o",
                str(output_dir),
                "--work-dir",
                str(work_dir),
            ],
        )
        assert result.exit_code == 0

        result = cli_runner.invoke(cli, ["run", str(config_file)])
        assert result.exit_code == 0

        output_files = list(output_dir.glob("*.nc"))
        assert len(output_files) > 0
        outputs.append(output_files[0])

    # Compare outputs
    ds1 = xr.open_dataset(outputs[0], engine="h5netcdf")
    ds2 = xr.open_dataset(outputs[1], engine="h5netcdf")

    # Check all data variables are identical
    for var in ds1.data_vars:
        if var not in ds2:
            continue
        np.testing.assert_array_equal(
            ds1[var].values,
            ds2[var].values,
            err_msg=f"{var} differs between runs",
        )

    ds1.close()
    ds2.close()


def test_golden_comparison(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
    cli_runner: CliRunner,
) -> None:
    """Compare current output against golden reference.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory.
    sample_disp_product : Path
        Mock DISP file.
    sample_static_los : Path
        Mock LOS file.
    sample_static_dem : Path
        Mock DEM file.
    sample_unr_data : tuple[Path, Path]
        (lookup_file, tenv8_dir).
    sample_algorithm_params : Path
        Sample algorithm parameters YAML file.
    cli_runner : CliRunner
        Click CLI test runner.

    """

    if cli is None:
        pytest.skip("Could not import CLI module")

    # Check for golden outputs
    golden_dir = Path(__file__).parent / "golden_output"
    if not golden_dir.exists():
        pytest.skip("Golden output directory not found")

    golden_files = list(golden_dir.glob("*.nc"))
    if not golden_files:
        pytest.skip("No golden output files found")

    # Run workflow
    lookup_file, tenv8_dir = sample_unr_data
    output_dir = tmp_path / "output"
    work_dir = tmp_path / "work"
    output_dir.mkdir()
    work_dir.mkdir()

    config_file = work_dir / "runconfig.yaml"
    result = cli_runner.invoke(
        cli,
        [
            "config",
            "-d",
            str(sample_disp_product),
            "-ul",
            str(lookup_file),
            "-ud",
            str(tenv8_dir),
            "-uv",
            "0.2",
            "-ut",
            "variable",
            "--los-file",
            str(sample_static_los),
            "--dem-file",
            str(sample_static_dem),
            "-a",
            str(sample_algorithm_params),
            "-c",
            str(config_file),
            "--frame-id",
            "8882",
            "-o",
            str(output_dir),
            "--work-dir",
            str(work_dir),
        ],
    )
    assert result.exit_code == 0

    result = cli_runner.invoke(cli, ["run", str(config_file)])
    assert result.exit_code == 0

    output_files = list(output_dir.glob("*.nc"))
    assert len(output_files) > 0

    # Compare against golden
    _compare_outputs(output_files[0], golden_files[0])


def _validate_output(output_file: Path) -> None:
    """Validate output file structure and data.

    Parameters
    ----------
    output_file : Path
        Output NetCDF file to validate.

    """
    ds = xr.open_dataset(output_file, engine="h5netcdf")

    # Check dimensions exist
    assert "y" in ds.dims or "latitude" in ds.dims
    assert "x" in ds.dims or "longitude" in ds.dims

    # Check some data variables exist
    assert len(ds.data_vars) > 0, "No data variables in output"

    # Check for reasonable data
    for var in ds.data_vars:
        data = ds[var].values

        # Should not be all NaN
        assert not np.all(np.isnan(data)), f"{var} is all NaN"

        # Check finite values exist
        assert np.any(np.isfinite(data)), f"{var} has no finite values"

    # Check attributes
    assert len(ds.attrs) > 0, "No global attributes"

    ds.close()


def _compare_outputs(
    current_file: Path,
    golden_file: Path,
    rtol: float = 1e-5,
    atol: float = 1e-8,
) -> None:
    """Compare current output against golden reference.

    Parameters
    ----------
    current_file : Path
        Current output file.
    golden_file : Path
        Golden reference file.
    rtol : float
        Relative tolerance for comparison.
    atol : float
        Absolute tolerance for comparison.

    """
    ds_current = xr.open_dataset(current_file, engine="h5netcdf")
    ds_golden = xr.open_dataset(golden_file, engine="h5netcdf")

    # Compare each variable
    for var in ds_golden.data_vars:
        if var not in ds_current:
            continue

        current_data = ds_current[var].values
        golden_data = ds_golden[var].values

        # Check shapes match
        assert current_data.shape == golden_data.shape, f"{var} shape mismatch"

        # Compare values (handle NaNs)
        valid_mask = ~np.isnan(golden_data) & ~np.isnan(current_data)

        if not np.any(valid_mask):
            continue

        np.testing.assert_allclose(
            current_data[valid_mask],
            golden_data[valid_mask],
            rtol=rtol,
            atol=atol,
            err_msg=f"{var} values differ from golden reference",
        )

    ds_current.close()
    ds_golden.close()


@pytest.mark.integration
def test_different_grid_types(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
    cli_runner: CliRunner,
) -> None:
    """Test workflow with different UNR grid types.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory.
    sample_disp_product : Path
        Mock DISP file.
    sample_static_los : Path
        Mock LOS file.
    sample_static_dem : Path
        Mock DEM file.
    sample_unr_data : tuple[Path, Path]
        (lookup_file, tenv8_dir).
    sample_algorithm_params : Path
        Sample algorithm parameters YAML file.
    cli_runner : CliRunner
        Click CLI test runner.

    """
    if cli is None:
        pytest.skip("Could not import CLI module")

    lookup_file, tenv8_dir = sample_unr_data
    grid_types = ["constant", "variable"]

    for grid_type in grid_types:
        output_dir = tmp_path / grid_type / "output"
        work_dir = tmp_path / grid_type / "work"
        output_dir.mkdir(parents=True)
        work_dir.mkdir(parents=True)

        config_file = work_dir / "runconfig.yaml"

        result = cli_runner.invoke(
            cli,
            [
                "config",
                "-d",
                str(sample_disp_product),
                "-ul",
                str(lookup_file),
                "-ud",
                str(tenv8_dir),
                "-uv",
                "0.2",
                "-ut",
                grid_type,
                "--los-file",
                str(sample_static_los),
                "--dem-file",
                str(sample_static_dem),
                "-a",
                str(sample_algorithm_params),
                "-c",
                str(config_file),
                "--frame-id",
                "8882",
                "-o",
                str(output_dir),
                "--work-dir",
                str(work_dir),
            ],
        )
        assert result.exit_code == 0, f"Config failed for {grid_type}"

        result = cli_runner.invoke(cli, ["run", str(config_file)])
        assert result.exit_code == 0, f"Run failed for {grid_type}"

        output_files = list(output_dir.glob("*.nc"))
        assert len(output_files) > 0, f"No outputs for {grid_type}"


@pytest.mark.integration
def test_output_data_quality(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
    cli_runner: CliRunner,
) -> None:
    """Test output data meets quality requirements.

    Parameters
    ----------
    tmp_path : Path
        Temporary directory.
    sample_disp_product : Path
        Mock DISP file.
    sample_static_los : Path
        Mock LOS file.
    sample_static_dem : Path
        Mock DEM file.
    sample_unr_data : tuple[Path, Path]
        (lookup_file, tenv8_dir).
    sample_algorithm_params : Path
        Sample algorithm parameters YAML file.
    cli_runner : CliRunner
        Click CLI test runner.

    """
    if cli is None:
        pytest.skip("Could not import CLI module")

    lookup_file, tenv8_dir = sample_unr_data
    output_dir = tmp_path / "output"
    work_dir = tmp_path / "work"
    output_dir.mkdir()
    work_dir.mkdir()

    config_file = work_dir / "runconfig.yaml"
    result = cli_runner.invoke(
        cli,
        [
            "config",
            "-d",
            str(sample_disp_product),
            "-ul",
            str(lookup_file),
            "-ud",
            str(tenv8_dir),
            "-uv",
            "0.2",
            "-ut",
            "variable",
            "--los-file",
            str(sample_static_los),
            "--dem-file",
            str(sample_static_dem),
            "-a",
            str(sample_algorithm_params),
            "-c",
            str(config_file),
            "--frame-id",
            "8882",
            "-o",
            str(output_dir),
            "--work-dir",
            str(work_dir),
        ],
    )
    assert result.exit_code == 0

    result = cli_runner.invoke(cli, ["run", str(config_file)])
    assert result.exit_code == 0

    output_files = list(output_dir.glob("*.nc"))
    assert len(output_files) > 0

    # Quality checks
    ds = xr.open_dataset(output_files[0], engine="h5netcdf")

    for var in ds.data_vars:
        data = ds[var].values

        # NaN fraction should be reasonable
        nan_frac = np.isnan(data).sum() / data.size
        assert nan_frac < 0.9, f"{var} has {nan_frac:.1%} NaNs"

        # Finite values should exist
        finite_data = data[np.isfinite(data)]
        assert len(finite_data) > 0, f"{var} has no finite values"

        # Check for inf values
        assert not np.any(np.isinf(data)), f"{var} contains inf values"

    ds.close()
