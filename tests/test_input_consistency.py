"""The static layers and TROPO files must match the DISP grid and frame."""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

import h5py
import numpy as np
import pytest
import rasterio
import xarray as xr
from rasterio.crs import CRS
from rasterio.transform import from_origin

from cal_disp.prep.consistency import (
    check_input_geometry,
    check_tropo_coverage,
    disp_geometry,
    mask_los_nodata,
)
from cal_disp.product import DispProduct
from cal_disp.workflow import _load_los_bands, run_calibration

# Grid of sample_disp_product / sample_static_los / sample_static_dem
ORIGIN = (405000.0, 3778000.0)
SPACING = 100.0
SHAPE = (200, 200)
LOS_NAME = (
    "OPERA_L3_DISP-S1-STATIC_F{frame:05d}_20140403_S1A_v1.0_line_of_sight_enu.tif"
)
DEM_NAME = "OPERA_L3_DISP-S1-STATIC_F{frame:05d}_20140403_S1A_v1.0_dem.tif"


def _write_static(
    path: Path,
    layer: str = "LOS",
    shape: tuple[int, int] = SHAPE,
    origin: tuple[float, float] = ORIGIN,
    spacing: float = SPACING,
    epsg: int | None = 32611,
) -> Path:
    """Write a LOS (3 bands) or DEM (1 band) GeoTIFF on the given grid."""
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 3 if layer == "LOS" else 1
    data = np.full((count, *shape), 0.5, dtype=np.float32)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=shape[0],
        width=shape[1],
        count=count,
        dtype="float32",
        crs=CRS.from_epsg(epsg) if epsg else None,
        transform=from_origin(origin[0], origin[1], spacing, spacing),
    ) as dst:
        dst.write(data)
    return path


def _name(layer: str, frame: int = 8882) -> str:
    return (LOS_NAME if layer == "LOS" else DEM_NAME).format(frame=frame)


class TestCheckInputGeometry:
    def test_matching_inputs_pass(
        self, sample_disp_product, sample_static_los, sample_static_dem
    ):
        geometries = check_input_geometry(
            sample_disp_product, sample_static_los, sample_static_dem
        )

        assert list(geometries) == ["DISP", "LOS", "DEM"]
        for geom in geometries.values():
            assert geom.shape == SHAPE
            assert geom.frame_id == 8882
            assert geom.crs.to_epsg() == 32611
            assert (geom.transform.c, geom.transform.f) == ORIGIN
        assert "frame F08882, shape (200, 200), EPSG:32611" in (
            geometries["LOS"].describe()
        )

    def test_accepts_a_disp_product_and_no_dem(
        self, sample_disp_product, sample_static_los
    ):
        geometries = check_input_geometry(
            DispProduct.from_path(sample_disp_product), sample_static_los
        )
        assert list(geometries) == ["DISP", "LOS"]

    @pytest.mark.parametrize("layer", ["LOS", "DEM"])
    @pytest.mark.parametrize(
        "kwargs, expected",
        [
            ({"shape": (199, 200)}, r"shape \(199, 200\) != \(200, 200\)"),
            (
                {"origin": (ORIGIN[0] + SPACING, ORIGIN[1])},
                r"origin shifted by \(1, -?0\) pixels",
            ),
            (
                {"origin": (ORIGIN[0], ORIGIN[1] - 1.0)},
                r"origin shifted by \(0, 0\.01\) pixels",
            ),
            ({"spacing": 90.0, "shape": SHAPE}, r"pixel \(90\.0, -90\.0\)"),
            ({"epsg": 32612}, "CRS EPSG:32612 != EPSG:32611"),
            ({"epsg": None}, "CRS none != EPSG:32611"),
        ],
    )
    def test_mismatched_grid_raises(
        self, tmp_path, sample_disp_product, sample_static_los, layer, kwargs, expected
    ):
        bad = _write_static(tmp_path / "bad" / _name(layer), layer, **kwargs)
        los, dem = (bad, None) if layer == "LOS" else (sample_static_los, bad)

        with pytest.raises(ValueError, match=expected) as excinfo:
            check_input_geometry(sample_disp_product, los, dem)
        # The message names the offending file and the DISP product
        assert f"{layer} file {bad}" in str(excinfo.value)
        assert sample_disp_product.name in str(excinfo.value)

    @pytest.mark.parametrize("layer", ["LOS", "DEM"])
    def test_other_frame_with_the_same_grid_raises(
        self, tmp_path, sample_disp_product, sample_static_los, layer
    ):
        """Same shape, transform and CRS, but the layers of frame 8883."""
        other = _write_static(tmp_path / "other" / _name(layer, frame=8883), layer)
        los, dem = (other, None) if layer == "LOS" else (sample_static_los, other)

        with pytest.raises(ValueError, match="frame id 8883 != 8882"):
            check_input_geometry(sample_disp_product, los, dem)

    def test_all_differences_of_a_file_are_reported(
        self, tmp_path, sample_disp_product
    ):
        bad = _write_static(
            tmp_path / _name("LOS", frame=1), shape=(10, 10), epsg=4326, spacing=1.0
        )
        with pytest.raises(ValueError) as excinfo:
            check_input_geometry(sample_disp_product, bad)
        message = str(excinfo.value)
        for part in ("frame id 1 != 8882", "shape (10, 10)", "CRS EPSG:4326", "pixel"):
            assert part in message

    def test_difference_below_tolerance_passes(self, tmp_path, sample_disp_product):
        # 1e-5 m = 1e-7 pixel < 1e-6 pixel
        los = _write_static(
            tmp_path / _name("LOS"), origin=(ORIGIN[0] + 1e-5, ORIGIN[1])
        )
        check_input_geometry(sample_disp_product, los)
        with pytest.raises(ValueError, match="origin shifted"):
            check_input_geometry(sample_disp_product, los, pixel_tolerance=1e-8)

    def test_non_opera_name_skips_frame_check_with_warning(
        self, tmp_path, sample_disp_product, caplog
    ):
        los = _write_static(tmp_path / "smooth_los.tif")
        with caplog.at_level(logging.WARNING, logger="cal_disp.prep.consistency"):
            geometries = check_input_geometry(sample_disp_product, los)

        assert geometries["LOS"].frame_id is None
        assert "frame id not checked" in caplog.text

    def test_grid_from_coordinates_without_geotransform(
        self, tmp_path, sample_disp_product, sample_static_los
    ):
        """Without GeoTransform the grid comes from the pixel-centre x/y."""
        disp = tmp_path / "no_gt" / sample_disp_product.name
        disp.parent.mkdir()
        shutil.copy(sample_disp_product, disp)
        with h5py.File(disp, "a") as f:
            del f["spatial_ref"].attrs["GeoTransform"]

        # x/y are centres: the corner is half a pixel up-left of x[0], y[0]
        transform = disp_geometry(disp).transform
        assert (transform.c, transform.f) == (ORIGIN[0] - 50.0, ORIGIN[1] + 50.0)
        assert (transform.a, transform.e) == (SPACING, -SPACING)

        with pytest.raises(ValueError, match=r"origin shifted by \(0\.5, 0\.5\)"):
            check_input_geometry(disp, sample_static_los)
        los = _write_static(
            tmp_path / _name("LOS"), origin=(ORIGIN[0] - 50.0, ORIGIN[1] + 50.0)
        )
        check_input_geometry(disp, los)

    def test_disp_without_crs_raises(
        self, tmp_path, sample_disp_product, sample_static_los
    ):
        disp = tmp_path / "no_crs" / sample_disp_product.name
        disp.parent.mkdir()
        shutil.copy(sample_disp_product, disp)
        with h5py.File(disp, "a") as f:
            del f["spatial_ref"].attrs["crs_wkt"]

        with pytest.raises(ValueError, match="has no CRS"):
            check_input_geometry(disp, sample_static_los)


class TestTropoCoverage:
    @staticmethod
    def _tropo(path: Path, lat: tuple[float, float], lon: tuple[float, float]) -> Path:
        ds = xr.Dataset(
            {"wet_delay": (["latitude", "longitude"], np.zeros((5, 5)))},
            coords={
                "latitude": np.linspace(lat[1], lat[0], 5),
                "longitude": np.linspace(lon[0], lon[1], 5),
            },
        )
        ds.to_netcdf(path, engine="h5netcdf")
        return path

    def test_covering_file_passes(
        self, tmp_path, sample_disp_product, sample_static_los
    ):
        tropo = self._tropo(tmp_path / "tropo.nc", (33.0, 35.0), (-119.0, -117.0))
        check_tropo_coverage(sample_disp_product, [tropo])
        check_input_geometry(
            sample_disp_product, sample_static_los, tropo_files=[tropo]
        )

    def test_file_not_covering_the_frame_raises(
        self, tmp_path, sample_disp_product, sample_static_los
    ):
        # The frame spans about 33.96-34.14 N: this file starts at 34.05 N
        tropo = self._tropo(tmp_path / "tropo.nc", (34.05, 35.0), (-119.0, -117.0))
        with pytest.raises(ValueError, match="does not cover DISP file") as excinfo:
            check_input_geometry(
                sample_disp_product, sample_static_los, tropo_files=[tropo]
            )
        assert str(tropo) in str(excinfo.value)
        assert "S 34.050" in str(excinfo.value)

    def test_tropo_preparation_checks_coverage(
        self, tmp_path, sample_disp_product, sample_static_los, sample_static_dem
    ):
        from cal_disp.prep.tropo import prepare_troposphere_correction

        tropo = self._tropo(tmp_path / "tropo.nc", (40.0, 41.0), (-119.0, -117.0))
        with pytest.raises(ValueError, match="does not cover DISP file"):
            prepare_troposphere_correction(
                disp_file=sample_disp_product,
                dem_file=sample_static_dem,
                los_file=sample_static_los,
                reference_tropo_files=[tropo],
                secondary_tropo_files=[tropo],
                output_dir=tmp_path / "tropo_out",
            )


class TestLosNodata:
    def test_all_zero_vectors_become_nan(self):
        east = np.array([[0.0, -0.6], [0.0, -0.6]], dtype=np.float32)
        north = np.array([[0.0, -0.1], [0.1, -0.1]], dtype=np.float32)
        up = np.array([[0.0, 0.8], [0.0, 0.8]], dtype=np.float32)

        e, n, u = mask_los_nodata(east, north, up)

        for band in (e, n, u):
            assert band.dtype == np.float32
            assert np.isnan(band[0, 0])
            assert np.isfinite(band[0, 1]) and np.isfinite(band[1, 0])
        # A vector with one zero component is a valid look
        assert e[1, 0] == 0.0 and n[1, 0] == np.float32(0.1)
        # The inputs are not modified
        assert east[0, 0] == 0.0

    def test_nodata_value_in_any_band_becomes_nan(self):
        east = np.array([-9999.0, -0.6], dtype=np.float32)
        north = np.array([-0.1, -0.1], dtype=np.float32)
        up = np.array([0.8, 0.8], dtype=np.float32)

        e, n, u = mask_los_nodata(east, north, up, nodata=-9999.0)

        assert np.isnan([e[0], n[0], u[0]]).all()
        assert np.isfinite([e[1], n[1], u[1]]).all()

    def test_valid_vectors_are_returned_unchanged(self):
        bands = tuple(np.full((3, 3), v, dtype=np.float32) for v in (-0.6, -0.1, 0.8))
        out = mask_los_nodata(*bands, nodata=0.0)
        assert all(o is b for o, b in zip(out, bands))

    def test_load_los_bands_masks_nodata(self, tmp_path):
        """The workflow's LOS reader returns NaN, not 0, outside the swath."""
        path = tmp_path / _name("LOS")
        data = np.full((3, 20, 30), 0.5, dtype=np.float32)
        data[:, :5, :] = 0.0  # outside the swath
        with rasterio.open(
            path,
            "w",
            driver="GTiff",
            height=20,
            width=30,
            count=3,
            dtype="float32",
            crs=CRS.from_epsg(32611),
            transform=from_origin(*ORIGIN, SPACING, SPACING),
            nodata=0.0,
        ) as dst:
            dst.write(data)

        east, north, up = _load_los_bands(path)

        for band in (east, north, up):
            assert np.isnan(band[:5]).all()
            np.testing.assert_array_equal(band[5:], 0.5)


class TestWorkflowRejectsInconsistentInputs:
    """run_calibration stops before any processing."""

    def _run(self, tmp_path, disp, unr, los, dem=None):
        lookup_file, tenv8_dir = unr
        return run_calibration(
            disp_file=disp,
            unr_grid_latlon_file=lookup_file,
            unr_timeseries_dir=tenv8_dir,
            output_dir=tmp_path / "out",
            los_file=los,
            dem_file=dem,
        )

    def test_shifted_los(self, tmp_path, sample_disp_product, sample_unr_data):
        los = _write_static(
            tmp_path / "shifted" / _name("LOS"), origin=(ORIGIN[0] + 300.0, ORIGIN[1])
        )
        with pytest.raises(ValueError, match="LOS file .* is not on the grid of DISP"):
            self._run(tmp_path, sample_disp_product, sample_unr_data, los)
        # Nothing was computed or written
        assert not (tmp_path / "out" / "scratch" / "gnss").exists()
        assert not list((tmp_path / "out").glob("*.nc"))

    def test_los_of_another_frame(self, tmp_path, sample_disp_product, sample_unr_data):
        los = _write_static(tmp_path / "other" / _name("LOS", frame=8883))
        with pytest.raises(ValueError, match="frame id 8883 != 8882"):
            self._run(tmp_path, sample_disp_product, sample_unr_data, los)

    def test_dem_in_another_crs(
        self, tmp_path, sample_disp_product, sample_unr_data, sample_static_los
    ):
        dem = _write_static(tmp_path / "dem" / _name("DEM"), "DEM", epsg=4326)
        with pytest.raises(ValueError, match="DEM file .* CRS EPSG:4326 != EPSG:32611"):
            self._run(
                tmp_path, sample_disp_product, sample_unr_data, sample_static_los, dem
            )
