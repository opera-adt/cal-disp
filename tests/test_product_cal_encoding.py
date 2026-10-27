"""DISP-CAL product writer: chunking and compression of the raster layers."""

from __future__ import annotations

import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest
import xarray as xr
from cal_product_helpers import write_product

from cal_disp.product import CalProduct
from cal_disp.product.output._utils import build_encoding


def test_build_encoding_rasters_only():
    ds = xr.Dataset(
        {
            "cal3d": (["time", "y", "x"], np.zeros((1, 300, 400), dtype=np.float32)),
            "cal2d": (["y", "x"], np.zeros((300, 400), dtype=np.float64)),
            "flag": (["y", "x"], np.zeros((300, 400), dtype=np.uint8)),
            "name": ((), "scalar string"),
            "spatial_ref": ((), 0),
        },
        coords={"time": [0], "y": np.arange(300.0), "x": np.arange(400.0)},
    )
    enc = build_encoding(ds)

    assert enc["cal3d"]["chunksizes"] == (1, 256, 256)
    assert enc["cal2d"]["chunksizes"] == (256, 256)
    assert enc["cal3d"]["zlib"] and enc["cal3d"]["shuffle"]
    assert enc["cal3d"]["complevel"] == 4
    assert enc["cal2d"]["dtype"] == "float32"
    assert np.isnan(enc["cal2d"]["_FillValue"])
    assert "_FillValue" not in enc["flag"]  # integer layer keeps its own
    assert "name" not in enc and "spatial_ref" not in enc
    assert enc["x"] == enc["y"] == enc["time"] == {"_FillValue": None}
    plain = build_encoding(ds, compression=False)
    assert not plain["cal3d"]["zlib"] and not plain["cal3d"]["shuffle"]
    assert plain["cal3d"]["chunksizes"] == (1, 256, 256)


@pytest.mark.parametrize("compression", [True, False])
def test_product_rasters_chunked_and_compressed(
    tmp_path: Path, sample_disp_product: Path, compression: bool
):
    product, cal = write_product(
        sample_disp_product, tmp_path / "out", compression=compression
    )

    with h5py.File(product.path) as f:
        for name in ("calibration", "calibration_std"):
            d = f[name]
            assert d.dtype == np.float32
            assert d.chunks == (1, 200, 200)  # (1, 256, 256) capped by the grid
            assert (d.compression == "gzip") is compression
            assert (d.compression_opts == 4) is compression
            assert d.shuffle is compression
            assert np.isnan(d.attrs["_FillValue"])
        aux = f["auxiliary/up_down"]
        assert aux.chunks == (1, 4, 4)
        assert (aux.compression == "gzip") is compression
        # coordinates, scalars and strings: no filters
        for name in ("x", "y", "time", "identification/frame_id"):
            assert f[name].compression is None
            assert not f[name].shuffle

    # Data round-trips bit-identically (NaNs included)
    with xr.open_dataset(product.path, engine="h5netcdf") as ds:
        written = ds["calibration"].values
        assert written.dtype == np.float32
        assert written.tobytes() == cal.tobytes()
        std = ds["calibration_std"].values
        assert std.tobytes() == (np.abs(cal) / 10).tobytes()


def test_compression_makes_file_smaller(tmp_path: Path, sample_disp_product: Path):
    zipped, _ = write_product(sample_disp_product, tmp_path / "z", compression=True)
    plain, _ = write_product(sample_disp_product, tmp_path / "p", compression=False)
    assert zipped.path.stat().st_size < 0.95 * plain.path.stat().st_size


def test_compression_plumbed_from_runconfig(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
):
    """output_options.compression reaches the writer through the PGE runconfig."""
    from cal_disp.cli.config import create_config
    from cal_disp.cli.run import run_main

    lookup_file, tenv8_dir = sample_unr_data
    outputs = {}
    for compression in (True, False):
        name = "zip" if compression else "plain"
        config_file = create_config(
            disp_file=sample_disp_product,
            frame_id=8882,
            unr_grid_latlon_file=lookup_file,
            unr_timeseries_dir=tenv8_dir,
            unr_grid_version="0.2",
            unr_grid_type="constant",
            algorithm_params_file=sample_algorithm_params,
            los_file=sample_static_los,
            dem_file=sample_static_dem,
            output_dir=tmp_path / name / "out",
            work_dir=tmp_path / name / "work",
            compression=compression,
        )
        outputs[compression] = run_main(config_file)

    with h5py.File(outputs[True]) as f:
        assert f["calibration"].compression == "gzip"
        assert f["calibration"].chunks == (1, 200, 200)
    with h5py.File(outputs[False]) as f:
        assert f["calibration"].compression is None
        assert f["calibration"].chunks == (1, 200, 200)


def test_from_path_defaults_to_compression(tmp_path: Path, sample_disp_product: Path):
    product, _ = write_product(sample_disp_product, tmp_path / "out")
    copy = tmp_path / "copy" / product.filename
    copy.parent.mkdir()
    shutil.copy(product.path, copy)
    assert CalProduct.from_path(copy).compression is True
