"""Depth of ``cal-disp validate``: every kind of product difference is caught."""

from __future__ import annotations

import logging
import shutil
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pytest
import xarray as xr
from click.testing import CliRunner

from cal_disp.cli.validate import validate_cli
from cal_disp.product import CalProduct, DispProduct
from cal_disp.validate import DEFAULT_TOLERANCE, compare_cal_products


def _build_product(
    disp_file: Path,
    out_dir: Path,
    cal_dtype: type = np.float32,
    with_auxiliary: bool = True,
) -> Path:
    """Write a small but complete DISP-CAL product (all four groups)."""
    disp = DispProduct.from_path(disp_file)
    with disp.open_dataset() as ds:
        ny, nx = ds["displacement"].shape
        coords = {"time": ds.time.values, "y": ds.y.values, "x": ds.x.values}
        spatial_ref = ds["spatial_ref"].copy()

    rng = np.random.default_rng(0)
    cal = (rng.normal(size=(ny, nx)) * 1e-3).astype(np.float32)
    cal[:10] = np.nan  # a NaN border, as in a real product
    std = np.abs(rng.normal(size=(ny, nx)) * 1e-4).astype(np.float32)
    std[:, :20] = 0.0  # clipped zeros, as in a real product

    def _da(data: np.ndarray, name: str) -> xr.DataArray:
        return xr.DataArray(
            data[np.newaxis].astype(cal_dtype),
            coords=coords,
            dims=["time", "y", "x"],
            attrs={"units": "meters", "long_name": name},
        )

    product = CalProduct.create(
        calibration=_da(cal, "calibration_correction"),
        disp_product=disp,
        output_dir=out_dir,
        calibration_std=_da(std, "calibration_uncertainty"),
        spatial_ref=spatial_ref,
        global_metadata={"gnss_reference_epoch": "2022.0274"},
        version="1.0",
    )
    product.add_identification(
        calibration_reference_name="UNR gridded data",
        calibration_reference_version="0.3",
        calibration_reference_type="constant",
        calibration_reference_reference_frame="IGS20",
        source_data_file_list=[disp_file.name],
        source_calibration_file_list=["grid_latlon_lookup_v0.3.txt"],
        source_data_access="https://datapool.asf.alaska.edu/",
        source_data_dem_name="Copernicus DEM GLO-30",
        source_data_satellite_names=["Sentinel-1A"],
        source_data_imaging_geometry="right_looking",
        source_data_x_spacing=100.0,
        source_data_y_spacing=100.0,
        static_layers_data_access="https://example.com/static_layers",
        absolute_orbit_number=41457,
        track_number=34,
        instrument_name="C-SAR",
        look_direction="right",
        radar_band="C",
        orbit_pass_direction="ascending",
        bounding_polygon="POLYGON((0 0, 1 0, 1 1, 0 1, 0 0))",
        product_bounding_box="(0, 0, 1, 1)",
        product_sample_spacing="100.0m",
        product_data_access="https://example.com/products",
        processing_facility="JPL",
        nodata_pixel_count=int(np.isnan(cal).sum()),
        ceos_number_of_input_granules=1,
        processing_start_datetime=datetime.now(tz=timezone.utc),
    )
    product.add_metadata(
        algorithm_parameters_yaml="calibration_options:\n  grid_type: constant\n",
        platform_id="S1A",
        source_data_software_disp_version="1.0",
        cal_disp_software_version="0.2",
        venti_software_version="0.1",
        pge_runconfig=f"disp_file: {disp_file}\n",
    )
    if with_auxiliary:
        coarse = {
            "time": coords["time"],
            "y": coords["y"][::50],
            "x": coords["x"][::50],
        }
        shape = (1, len(coarse["y"]), len(coarse["x"]))
        product.add_auxiliary(
            model_3d={
                name: xr.DataArray(
                    np.zeros(shape, dtype=np.float32),
                    coords=coarse,
                    dims=["time", "y", "x"],
                    attrs={"units": "meters/year", "long_name": name},
                )
                for name in ("north_south", "east_west", "up_down")
            },
            spatial_ref=spatial_ref,
        )
    return product.path


@pytest.fixture
def reference_product(sample_disp_product: Path, tmp_path: Path) -> Path:
    return _build_product(sample_disp_product, tmp_path / "reference")


@pytest.fixture
def test_product(reference_product: Path, tmp_path: Path) -> Path:
    """A byte-identical copy of the reference, to be mutated by each test."""
    out = tmp_path / "test" / reference_product.name
    out.parent.mkdir()
    shutil.copy(reference_product, out)
    return out


def _validate(ref: Path, test: Path, caplog, **kwargs) -> tuple[bool, str]:
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="cal_disp.validate"):
        ok = compare_cal_products(ref, test, **kwargs)
    return ok, caplog.text


def test_identical_copy_passes(reference_product, test_product, caplog):
    ok, log = _validate(reference_product, test_product, caplog)
    assert ok, log
    # No browse image next to the test product is only a warning
    assert "no browse image" in log


def test_self_comparison_passes(reference_product, caplog):
    ok, log = _validate(reference_product, reference_product, caplog)
    assert ok, log


# One mutation at a time; each must be caught and named in the report


def test_catches_crs_change(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        wkt = f["spatial_ref"].attrs["crs_wkt"]
        f["spatial_ref"].attrs["crs_wkt"] = wkt.replace("zone 11N", "zone 12N").replace(
            '"central_meridian",-117', '"central_meridian",-111'
        )
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "spatial_ref: CRS EPSG:32612, expected EPSG:32611" in log


def test_catches_grid_shift(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["x"][:] = f["x"][:] + 1.0  # shift the grid by 1 m (1e-2 px)
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "x grid: origin 405001.0, step 100.0; expected origin 405000.0" in log


def test_catches_geotransform_mismatch(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["spatial_ref"].attrs["GeoTransform"] = "0 1 0 0 0 -1"
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "spatial_ref: GeoTransform '0 1 0 0 0 -1', expected '" in log

    # A reference without GeoTransform (products written before it was kept)
    with h5py.File(reference_product, "a") as f:
        del f["spatial_ref"].attrs["GeoTransform"]
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "spatial_ref: GeoTransform missing from the reference" in log

    with h5py.File(reference_product, "a") as f:
        f["spatial_ref"].attrs["GeoTransform"] = "1 1 0 0 0 -1"
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "GeoTransform '0 1 0 0 0 -1', expected '1 1 0 0 0 -1'" in log


def test_catches_units_attribute_change(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["calibration"].attrs["units"] = "millimeters"
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "calibration attribute 'units': 'millimeters', expected 'meters'" in log


def test_catches_global_attribute_change(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f.attrs["title"] = "Something else"
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "group attribute 'title'" in log


def test_catches_dtype_change(reference_product, test_product, caplog):
    # The writer always encodes rasters as float32, so change the on-disk
    # dtype directly: replace each raster by a float64 copy with its attrs.
    with h5py.File(test_product, "a") as f:
        for name in ("calibration", "calibration_std"):
            old = f[name]
            data, attrs = old[...].astype(np.float64), dict(old.attrs)
            scales = [[s[1].name for s in dim.items()] for dim in old.dims]
            del f[name]
            ds = f.create_dataset(name, data=data)
            ds.attrs.update(attrs)
            for dim, names in zip(ds.dims, scales):
                for scale_name in names:
                    dim.attach_scale(f[scale_name])
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "calibration: dtype float64, expected float32" in log
    assert "calibration_std: dtype float64, expected float32" in log


def test_catches_identification_value_change(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["identification/track_number"][()] = 0
        f["identification/product_version"][()] = "0.3."
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "track_number: 0, expected 34" in log
    assert "product_version: '0.3.', expected '1.0'" in log


def test_catches_metadata_value_change(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["metadata/platform_id"][()] = "S1B"
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "platform_id: 'S1B', expected 'S1A'" in log


def test_catches_missing_group(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        del f["auxiliary"]
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "[auxiliary]" in log
    assert "group missing from the test product" in log


def test_catches_missing_identification_group(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        del f["identification"]
    ok, log = _validate(reference_product, test_product, caplog, group="main")
    assert not ok
    assert "[identification]" in log


def test_catches_missing_variable(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        del f["calibration_std"]
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "missing from test: ['calibration_std']" in log


def test_catches_value_beyond_tolerance(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["calibration"][0, 50, 50] += 1e-4
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "calibration: 1 of" in log
    assert "first at index (0, 50, 50)" in log


def test_value_within_tolerance_passes(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["calibration"][0, 50, 50] += 1e-8
    ok, log = _validate(reference_product, test_product, caplog)
    assert ok, log


def test_reference_zeros_are_compared(reference_product, test_product, caplog):
    """Pixels where the reference is exactly 0 are compared with atol."""
    with h5py.File(reference_product) as f:
        assert f["calibration_std"][0, 5, 5] == 0.0
    with h5py.File(test_product, "a") as f:
        f["calibration_std"][0, 5, 5] = 1e-3
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "calibration_std: 1 of" in log


def test_catches_nan_pattern_change(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["calibration"][0, 100, 100] = np.nan
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "calibration: NaN pattern differs" in log


def test_rtol_and_atol_are_separate(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["calibration"][0, 50, 50] += 1e-4
    ok, _ = _validate(reference_product, test_product, caplog, rtol=0, atol=1e-3)
    assert ok
    ok, _ = _validate(reference_product, test_product, caplog, rtol=1e-3, atol=0)
    assert not ok  # relative change is ~1e-1
    ok, _ = _validate(reference_product, test_product, caplog, tolerance=1e-3)
    assert ok


def test_volatile_fields_are_ignored(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f.attrs["software_version"] = "9.9.9"
        f["identification/processing_start_datetime"][()] = "2030-01-01T00:00:00"
        f["metadata/cal_disp_software_version"][()] = "9.9.9"
        f["metadata/venti_software_version"][()] = "9.9.9"
        f["metadata/pge_runconfig"][()] = "disp_file: /somewhere/else.nc\n"
    ok, log = _validate(reference_product, test_product, caplog)
    assert ok, log


def test_catches_filename_metadata_change(reference_product, test_product, caplog):
    renamed = test_product.with_name(test_product.name.replace("_VV_", "_VH_"))
    test_product.rename(renamed)
    ok, log = _validate(reference_product, renamed, caplog)
    assert not ok
    assert "polarization: reference=VV, test=VH" in log


def test_group_selection(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["auxiliary/up_down"][0, 0, 0] = 1.0
    ok, _ = _validate(reference_product, test_product, caplog, group="main")
    assert ok
    ok, log = _validate(reference_product, test_product, caplog, group="auxiliary")
    assert not ok
    assert "up_down: 1 of" in log


# Browse image


def _write_png(path: Path, size: tuple[int, int]) -> None:
    from PIL import Image

    Image.new("RGB", size).save(path)


def test_browse_image_ok(reference_product, test_product, caplog):
    _write_png(test_product.with_suffix(".png"), (1024, 800))
    ok, log = _validate(reference_product, test_product, caplog)
    assert ok, log
    assert "no browse image" not in log


def test_browse_image_too_large_fails(reference_product, test_product, caplog):
    _write_png(test_product.with_suffix(".png"), (2049, 100))
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "is 2049x100 px; each side must be <= 2048 px" in log


def test_all_differences_reported_in_one_run(reference_product, test_product, caplog):
    with h5py.File(test_product, "a") as f:
        f["calibration"].attrs["units"] = "mm"
        f["identification/track_number"][()] = 0
        del f["auxiliary"]
    ok, log = _validate(reference_product, test_product, caplog)
    assert not ok
    assert "units" in log and "track_number" in log and "[auxiliary]" in log


# CLI


def test_cli_passes_and_fails_for_real_products(reference_product, test_product):
    runner = CliRunner()
    result = runner.invoke(
        validate_cli, [str(reference_product), str(test_product)], obj={"debug": False}
    )
    assert result.exit_code == 0, result.output

    with h5py.File(test_product, "a") as f:
        f["identification/frame_id"][()] = 1
    result = runner.invoke(
        validate_cli, [str(reference_product), str(test_product)], obj={"debug": False}
    )
    assert result.exit_code == 1


def test_cli_tolerance_options(reference_product, test_product):
    with h5py.File(test_product, "a") as f:
        f["calibration"][0, 50, 50] += 1e-4
    runner = CliRunner()
    args = [str(reference_product), str(test_product)]
    assert runner.invoke(validate_cli, args, obj={"debug": False}).exit_code == 1
    assert (
        runner.invoke(
            validate_cli, [*args, "--atol", "1e-3"], obj={"debug": False}
        ).exit_code
        == 0
    )
    assert (
        runner.invoke(
            validate_cli, [*args, "--tolerance", "1e-3"], obj={"debug": False}
        ).exit_code
        == 0
    )
    result = runner.invoke(
        validate_cli,
        [*args, "--tolerance", "1e-3", "--rtol", "1"],
        obj={"debug": False},
    )
    assert result.exit_code == 2
    assert "cannot be combined" in result.output


def test_one_default_tolerance():
    """The CLI default is the function default."""
    from cal_disp.cli.validate import validate_cli as cli

    assert DEFAULT_TOLERANCE == 1e-6
    for opt in cli.params:
        if opt.name in ("rtol", "atol", "tolerance"):
            assert opt.default is None  # resolved to DEFAULT_TOLERANCE at run time
