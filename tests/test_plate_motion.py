"""The plate_motion group and the GNSS displacement in /auxiliary."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest
import xarray as xr
from geepers.euler import plate_velocity_enu
from pyproj import Transformer

from cal_disp.config._algorithm import AlgorithmParameters
from cal_disp.validate import compare_cal_products
from cal_disp.workflow import AUXILIARY_STEP, _plate_motion_los, run_calibration

EPSG = 32611


def _grid(n: int = 400, step: float = 30.0):
    x = 405000.0 + step / 2 + step * np.arange(n)
    y = 3778000.0 - step / 2 - step * np.arange(n)
    return x, y


def test_plate_motion_matches_the_pole_at_every_pixel():
    """Coarse evaluation + interpolation reproduces a per-pixel ITRF2020-PMM
    prediction projected on the LOS and scaled by the span (to 1 micron)."""
    x, y = _grid()
    shape = (len(y), len(x))
    e = np.full(shape, -0.62, np.float32)
    n = np.full(shape, -0.11, np.float32)
    u = np.full(shape, 0.78, np.float32)
    e[:5, :5] = n[:5, :5] = 0  # no look
    dt = 0.5
    out = _plate_motion_los("NA", EPSG, (e, n, u), x, y, dt, block_rows=64)

    rows = np.array([10, 200, 399])
    cols = np.array([7, 250, 399])
    xx, yy = np.meshgrid(x[cols], y[rows])
    lon, lat = Transformer.from_crs(EPSG, 4326, always_xy=True).transform(xx, yy)
    ve, vn, _ = plate_velocity_enu(lon.ravel(), lat.ravel(), "NA")
    expect = (-0.62 * np.asarray(ve) - 0.11 * np.asarray(vn)) * dt / 1000.0
    got = out[np.ix_(rows, cols)].ravel()
    np.testing.assert_allclose(got, expect, atol=1e-6)
    assert np.isnan(out[:5, :5]).all()
    assert np.isfinite(out[5:, 5:]).all()


def test_plate_motion_scales_with_the_span_and_the_plate():
    x, y = _grid(50)
    los = tuple(np.full((50, 50), v, np.float32) for v in (-0.6, -0.1, 0.79))
    one = _plate_motion_los("NA", EPSG, los, x, y, 1.0)
    two = _plate_motion_los("NA", EPSG, los, x, y, 2.0)
    pacific = _plate_motion_los("PA", EPSG, los, x, y, 1.0)
    np.testing.assert_allclose(two, 2 * one, atol=1e-6)  # each rounded to 1 micron
    assert np.abs(pacific - one).min() > 1e-3  # NA and PA differ by cm/yr here


def _run(tmp_path, disp, los, unr, plate="NA"):
    lookup_file, tenv8_dir = unr
    params = AlgorithmParameters.create_default().with_overrides(
        {
            "calibration_options": {
                "apply_solid_earth_tide_correction": False,
                "downsample_factor": 2,
                "surface": {"cutoff_wavelength_meters": 5000.0},
                "frame": {"plate": plate},
            }
        }
    )
    return run_calibration(
        disp_file=disp,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / f"out_{plate}",
        los_file=los,
        algorithm_parameters=params,
    )


def test_product_carries_plate_motion_and_gnss_displacement(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    out = _run(tmp_path, sample_disp_product, sample_static_los, sample_unr_data)
    with xr.open_dataset(out) as main, xr.open_dataset(out, group="plate_motion") as pm:
        assert pm["plate_motion"].shape == main["calibration"].shape
        np.testing.assert_array_equal(pm.x.values, main.x.values)
        assert pm.attrs["plate"] == "NA"
        assert pm.attrs["plate_motion_model"] == "ITRF2020-PMM"
        span = pm.attrs["pair_span_years"]
        assert 0.5 < span < 0.6  # 2022-01-11 -> 2022-07-22
        assert pm["plate_motion"].attrs["units"] == "meters"
        assert pm["plate_motion"].attrs["grid_mapping"] == "spatial_ref"
        assert np.isfinite(pm["plate_motion"].values).mean() > 0.9
        assert "spatial_ref" in pm
    with xr.open_dataset(out, group="auxiliary") as aux:
        assert aux.attrs["pair_span_years"] == pytest.approx(span)
        assert "UNR grid" in aux.attrs["description"]
        for comp in ("east_west", "north_south", "up_down"):
            assert np.isfinite(aux[comp].values).all()
            assert f"{comp}_std" in aux
        assert np.abs(aux["east_west"].values).max() > 0  # no longer a placeholder
        assert aux["east_west"].shape[1:] == (
            len(aux.y),
            len(aux.x),
        )
    with xr.open_dataset(sample_disp_product) as ds:
        assert len(aux.x) == len(ds.x.values[::AUXILIARY_STEP])


def test_frame_plate_reaches_the_product(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    na = _run(tmp_path, sample_disp_product, sample_static_los, sample_unr_data, "NA")
    pa = _run(tmp_path, sample_disp_product, sample_static_los, sample_unr_data, "PA")
    with xr.open_dataset(na) as a, xr.open_dataset(pa) as b:
        # calibration does not depend on the plate: plate_motion is not part of it
        np.testing.assert_array_equal(a["calibration"].values, b["calibration"].values)
    with xr.open_dataset(pa, group="plate_motion") as pm:
        assert pm.attrs["plate"] == "PA"


def _validate(ref: Path, test: Path, caplog) -> tuple[bool, str]:
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="cal_disp.validate"):
        ok = compare_cal_products(ref, test)
    return ok, caplog.text


def test_validate_compares_the_plate_motion_group(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    caplog,
):
    ref = _run(tmp_path, sample_disp_product, sample_static_los, sample_unr_data)
    test = tmp_path / "test" / ref.name
    test.parent.mkdir()
    shutil.copy(ref, test)
    for png in ref.parent.glob("*.png"):
        shutil.copy(png, test.parent / png.name)
    ok, _ = _validate(ref, test, caplog)
    assert ok

    with h5py.File(test, "a") as f:
        f["plate_motion/plate_motion"][0, 50, 50] += 0.01
    ok, log = _validate(ref, test, caplog)
    assert not ok
    assert "[plate_motion]" in log

    with h5py.File(test, "a") as f:
        del f["plate_motion"]
    ok, log = _validate(ref, test, caplog)
    assert not ok
    assert "group missing from the test product" in log


def test_plate_motion_from_the_los_file_matches_in_memory(sample_static_los: Path):
    """The run reads the LOS file in row windows after the fit; same values as
    the in-memory bands, with the file's no-look pixels NaN."""
    import rasterio

    with rasterio.open(sample_static_los) as ds:
        bands = ds.read((1, 2, 3)).astype(np.float32)
        t, epsg = ds.transform, ds.crs.to_epsg()
        x = t.c + t.a * (np.arange(ds.width) + 0.5)
        y = t.f + t.e * (np.arange(ds.height) + 0.5)
    from_file = _plate_motion_los(
        "NA", epsg, sample_static_los, x, y, 0.7, block_rows=37
    )
    from_mem = _plate_motion_los("NA", epsg, tuple(bands), x, y, 0.7, block_rows=37)
    np.testing.assert_array_equal(np.isnan(from_file), np.isnan(from_mem))
    np.testing.assert_allclose(from_file, from_mem, atol=0, equal_nan=True)
    assert np.isfinite(from_file).mean() > 0.9
