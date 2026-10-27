"""Runconfig `fit_mask_file`: a stable-pixel fit mask not packaged in Venti."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import rasterio
import xarray as xr
import yaml  # type: ignore[import-untyped]
from rasterio.transform import from_origin

from cal_disp.cli.config import config_cli
from cal_disp.config._algorithm import AlgorithmParameters
from cal_disp.workflow import run_calibration


@pytest.fixture
def fit_mask_file(tmp_path: Path, sample_disp_product: Path) -> Path:
    """uint8 mask on the DISP grid that leaves the left third out of the fit."""
    with xr.open_dataset(sample_disp_product) as ds:
        x, y = ds.x.values, ds.y.values
        crs = ds["spatial_ref"].attrs["crs_wkt"]
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    transform = from_origin(x[0] - dx / 2, y[0] - dy / 2, dx, -dy)
    mask = np.ones((len(y), len(x)), dtype=np.uint8)
    mask[:, : len(x) // 3] = 0
    path = tmp_path / "F08882_test_stable.tif"
    with rasterio.open(
        path, "w", driver="GTiff", height=len(y), width=len(x), count=1,
        dtype="uint8", crs=crs, transform=transform,
    ) as dst:  # fmt: skip
        dst.write(mask, 1)
    return path


def _params(**surface) -> AlgorithmParameters:
    return AlgorithmParameters.create_default().with_overrides(
        {
            "calibration_options": {
                "apply_solid_earth_tide_correction": False,
                "downsample_factor": 2,
                "surface": {"cutoff_wavelength_meters": 5000.0, **surface},
            }
        }
    )


def _run(tmp_path, name, disp, los, unr, params, **kw) -> np.ndarray:
    lookup_file, tenv8_dir = unr
    out = run_calibration(
        disp_file=disp,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / name,
        los_file=los,
        algorithm_parameters=params,
        **kw,
    )
    with xr.open_dataset(out) as ds:
        return ds["calibration"].values


def test_fit_mask_file_equals_the_same_mask_packaged(
    tmp_path, sample_disp_product, sample_static_los, sample_unr_data, fit_mask_file,
    monkeypatch,
):  # fmt: skip
    """A mask from the runconfig fits exactly like the same mask by packaged name,
    and differs from no mask at all."""
    import venti.frames

    args = (sample_disp_product, sample_static_los, sample_unr_data)
    from_file = _run(tmp_path, "file", *args, _params(), fit_mask_file=fit_mask_file)
    monkeypatch.setattr(venti.frames, "FIT_MASK_DIR", fit_mask_file.parent)
    packaged = _run(tmp_path, "pkg", *args, _params(fit_mask=fit_mask_file.name))
    without = _run(tmp_path, "none", *args, _params())
    np.testing.assert_array_equal(from_file, packaged)
    assert np.nanmax(np.abs(from_file - without)) > 0


def test_fit_mask_file_and_packaged_name_are_exclusive(
    tmp_path, sample_disp_product, sample_static_los, sample_unr_data, fit_mask_file
):
    with pytest.raises(ValueError, match="both set"):
        _run(
            tmp_path, "both", sample_disp_product, sample_static_los, sample_unr_data,
            _params(fit_mask="F08882_stable_coh07.tif"), fit_mask_file=fit_mask_file,
        )  # fmt: skip


def test_fit_mask_file_needs_loclin(
    tmp_path, sample_disp_product, sample_static_los, sample_unr_data, fit_mask_file
):
    with pytest.raises(ValueError, match="loclin"):
        _run(
            tmp_path, "gamma", sample_disp_product, sample_static_los, sample_unr_data,
            _params(method="windowed_plane"), fit_mask_file=fit_mask_file,
        )  # fmt: skip


def test_config_cli_writes_fit_mask_file(
    cli_runner, tmp_path, sample_disp_product, sample_unr_grid_latlon,
    sample_unr_timeseries_dir, sample_algorithm_params, sample_static_layers,
    fit_mask_file,
):  # fmt: skip
    los_file, dem_file = sample_static_layers
    config_file = tmp_path / "runconfig.yaml"
    result = cli_runner.invoke(
        config_cli,
        [
            "--config-file", str(config_file),
            "--disp-file", str(sample_disp_product),
            "--unr-grid-latlon", str(sample_unr_grid_latlon),
            "--unr-grid-dir", str(sample_unr_timeseries_dir),
            "--frame-id", "8882",
            "--unr-grid-version", "0.2",
            "--algorithm-params", str(sample_algorithm_params),
            "--los-file", str(los_file),
            "--dem-file", str(dem_file),
            "--work-dir", str(tmp_path),
            "--output-dir", str(tmp_path / "out"),
            "--fit-mask-file", str(fit_mask_file),
        ],
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    rc = yaml.safe_load(config_file.read_text())["cal_disp_workflow"]
    assert rc["static_ancillary_group"]["fit_mask_file"] == str(fit_mask_file.resolve())
