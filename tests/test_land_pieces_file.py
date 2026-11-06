"""Runconfig `land_pieces_file`: breaklines (Venti surface.regions = 'pieces')."""

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


def _write_labels(path: Path, product: Path, dtype: str, split: bool) -> Path:
    with xr.open_dataset(product) as ds:
        x, y = ds.x.values, ds.y.values
        crs = ds["spatial_ref"].attrs["crs_wkt"]
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    transform = from_origin(x[0] - dx / 2, y[0] - dy / 2, dx, -dy)
    labels = np.ones((len(y), len(x)), dtype=dtype)
    if split:
        labels[:, : len(x) // 3] = 2  # the left third is a piece of its own
    with rasterio.open(
        path, "w", driver="GTiff", height=len(y), width=len(x), count=1,
        dtype=dtype, crs=crs, transform=transform,
    ) as dst:  # fmt: skip
        dst.write(labels, 1)
    return path


@pytest.fixture
def land_pieces_file(tmp_path: Path, sample_disp_product: Path) -> Path:
    return _write_labels(
        tmp_path / "F08882_land_pieces.tif", sample_disp_product, "uint16", True
    )


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


def test_pieces_give_the_split_piece_its_own_surface(
    tmp_path, sample_disp_product, sample_static_los, sample_unr_data, land_pieces_file
):
    args = (sample_disp_product, sample_static_los, sample_unr_data)
    params = _params(regions="pieces", piece_min_pixels=10)
    pieces = _run(tmp_path, "pieces", *args, params, land_pieces_file=land_pieces_file)
    none = _run(tmp_path, "none", *args, _params(), land_pieces_file=land_pieces_file)
    assert np.nanmax(np.abs(pieces - none)) > 0


def test_one_piece_is_the_frame_surface(
    tmp_path, sample_disp_product, sample_static_los, sample_unr_data
):
    one = _write_labels(tmp_path / "one.tif", sample_disp_product, "uint8", False)
    args = (sample_disp_product, sample_static_los, sample_unr_data)
    pieces = _run(
        tmp_path, "one", *args, _params(regions="pieces"), land_pieces_file=one
    )
    none = _run(tmp_path, "none", *args, _params())
    np.testing.assert_allclose(pieces, none, atol=1e-7, equal_nan=True)


def test_pieces_need_the_file(
    tmp_path, sample_disp_product, sample_static_los, sample_unr_data
):
    with pytest.raises(ValueError, match="land_pieces_file"):
        _run(
            tmp_path, "nofile", sample_disp_product, sample_static_los, sample_unr_data,
            _params(regions="pieces"),
        )  # fmt: skip


def test_labels_must_be_integers_on_the_grid(
    tmp_path, sample_disp_product, sample_static_los, sample_unr_data
):
    bad = _write_labels(tmp_path / "bad.tif", sample_disp_product, "float32", True)
    with pytest.raises(ValueError, match="dtype"):
        _run(
            tmp_path, "bad", sample_disp_product, sample_static_los, sample_unr_data,
            _params(regions="pieces"), land_pieces_file=bad,
        )  # fmt: skip


def test_config_cli_writes_land_pieces_file(
    cli_runner, tmp_path, sample_disp_product, sample_unr_grid_latlon,
    sample_unr_timeseries_dir, sample_algorithm_params, sample_static_layers,
    land_pieces_file,
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
            "--land-pieces-file", str(land_pieces_file),
        ],
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    rc = yaml.safe_load(config_file.read_text())["cal_disp_workflow"]
    group = rc["static_ancillary_group"]
    assert group["land_pieces_file"] == str(land_pieces_file.resolve())
