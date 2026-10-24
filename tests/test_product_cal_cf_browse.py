"""DISP-CAL product: CF/attribute hygiene and the browse image."""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import xarray as xr
from cal_product_helpers import write_product
from PIL import Image

from cal_disp._version import __version__
from cal_disp.browse_image import (
    _resize_to_max_pixel_dim,
    make_browse_image_from_arr,
    make_browse_image_from_nc,
    symmetric_limits,
)
from cal_disp.product import CalProduct


def test_root_attrs(tmp_path: Path, sample_disp_product: Path):
    product, _ = write_product(sample_disp_product, tmp_path / "out")

    with xr.open_dataset(product.path, engine="h5netcdf") as ds:
        attrs = ds.attrs
    assert attrs["mission_name"] == "OPERA"
    assert "mision_name" not in attrs
    assert attrs["spatial_resolution"] == "100 meters"
    assert attrs["temporal_resolution"] == "192 days"  # 2022-01-11 -> 2022-07-22
    assert "TBD" not in attrs.values()
    assert attrs["history"].endswith(f"created by cal_disp {__version__}")
    # the software version is stored once, in /metadata
    assert "software_version" not in attrs
    meta = CalProduct.from_path(product.path).open_metadata()
    assert "cal_disp_software_version" in meta


def test_coordinate_variables_follow_cf(tmp_path: Path, sample_disp_product: Path):
    product, _ = write_product(sample_disp_product, tmp_path / "out")

    for group in (None, "auxiliary"):
        with xr.open_dataset(product.path, group=group, engine="h5netcdf") as ds:
            for name, standard_name in (
                ("x", "projection_x_coordinate"),
                ("y", "projection_y_coordinate"),
            ):
                assert ds[name].attrs["standard_name"] == standard_name
                assert ds[name].attrs["units"] == "m"
                assert ds[name].attrs["long_name"]
            assert ds["time"].attrs["standard_name"] == "time"
            for var in ds.data_vars:
                if var == "spatial_ref":
                    continue
                assert ds[var].attrs["grid_mapping"] == "spatial_ref"
                assert ds[var].attrs["units"] == "meters"
                assert ds[var].dims == ("time", "y", "x")

    with h5py.File(product.path) as f:
        # CF: coordinate variables must not have _FillValue
        for path in ("x", "y", "time", "auxiliary/x", "auxiliary/y", "auxiliary/time"):
            assert "_FillValue" not in f[path].attrs, path
        # no hand-set "coordinates" attribute (was "y x" on (time, y, x) layers)
        for path in ("calibration", "calibration_std", "auxiliary/up_down"):
            assert "coordinates" not in f[path].attrs, path
            assert f[path].attrs["grid_mapping"] == "spatial_ref"
            assert "spatial_ref" in f[path].parent


def test_resize_uses_min_ratio_and_keeps_input():
    arr = np.random.default_rng(0).normal(size=(300, 500)).astype(np.float32)
    arr[:, :20] = np.nan
    before = arr.copy()

    out = _resize_to_max_pixel_dim(arr, 100)

    assert arr.tobytes() == before.tobytes()  # not mutated
    # the *largest* dimension is scaled to the limit (was the smallest)
    assert out.shape == (60, 100)
    assert np.isnan(out[:, :3]).all()  # NaN band preserved
    assert np.isfinite(out[:, 5:]).all()
    # no zero-bleeding at the NaN edge: same spread as the interior
    assert out[:, 4:6].std() > 0.5 * out[:, 50:].std()
    # small arrays are never upsampled
    assert _resize_to_max_pixel_dim(arr, 2048).shape == arr.shape


def test_resize_full_frame_shape_within_limit():
    """The F08882 frame (7733 x 9464) gave a 2506 x 2048 PNG before the fix."""
    out = _resize_to_max_pixel_dim(np.zeros((7733, 9464), dtype=np.float32), 2048)
    assert max(out.shape) <= 2048
    assert out.shape == (1673, 2048)


def test_symmetric_limits():
    arr = np.concatenate([np.linspace(-0.083, 0.155, 1000), [np.nan, 5.0]])
    vmin, vmax = symmetric_limits(arr, 98)
    assert vmin == -vmax
    assert 0.14 < vmax < 0.156  # outlier ignored; span not clipped at 0.10
    assert symmetric_limits(np.full(4, np.nan)) == (-1e-3, 1e-3)


def test_browse_image_size_and_limits(tmp_path: Path, sample_disp_product: Path):
    product, cal = write_product(sample_disp_product, tmp_path / "out")
    png = tmp_path / "browse.png"

    vmin, vmax = make_browse_image_from_nc(png, product.path, max_dim_allowed=128)

    with Image.open(png) as im:
        assert im.size == (128, 128)
        assert im.mode == "RGBA"
        alpha = np.asarray(im)[..., 3]
    assert (alpha[:5, :5] == 0).all()  # NaN block transparent
    assert (alpha[20:, 20:] == 255).all()
    assert vmin == -vmax and 0.09 < vmax < 0.2

    # explicit limits and mask honoured; input array untouched
    arr = cal[0].copy()
    mask = np.ones(arr.shape, dtype=bool)
    mask[:, :100] = False
    lim = make_browse_image_from_arr(png, arr, mask, 64, vmin=-0.1, vmax=0.1)
    assert lim == (-0.1, 0.1)
    assert arr.tobytes() == cal[0].tobytes()
    with Image.open(png) as im:
        assert im.size == (64, 64)
        assert (np.asarray(im)[..., 3][:, :30] == 0).all()
