"""Tests for TropoProduct class."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

from cal_disp.product import TropoProduct


class TestTropoProductParsing:
    """Tests for filename parsing."""

    def test_from_path(self):
        """Should parse filename correctly."""
        filename = (
            "OPERA_L4_TROPO-ZENITH_20220111T000000Z_20220111T120000Z_HRRR_v1.0.nc"
        )
        product = TropoProduct.from_path(filename)

        assert product.date == datetime(2022, 1, 11, 0, 0, 0)
        assert product.production_date == datetime(2022, 1, 11, 12, 0, 0)
        assert product.model == "HRRR"
        assert product.version == "1.0"

    def test_from_path_invalid(self):
        """Should reject invalid filename."""
        with pytest.raises(ValueError, match="does not match"):
            TropoProduct.from_path("invalid_file.nc")

    def test_validates_production_date(self):
        """Should reject production before model date."""
        with pytest.raises(ValueError, match="cannot be before"):
            TropoProduct(
                path=Path("test.nc"),
                date=datetime(2022, 1, 11),
                production_date=datetime(2022, 1, 10),
                model="HRRR",
                version="1.0",
            )


class TestTropoProductDataAccess:
    """Tests for data access methods."""

    def test_open_dataset(self, sample_tropo_product: Path):
        """Should open dataset."""
        product = TropoProduct.from_path(sample_tropo_product)
        ds = product.open_dataset()

        assert "wet_delay" in ds
        assert "hydrostatic_delay" in ds
        assert "height" in ds.dims
        assert "latitude" in ds.dims
        assert "longitude" in ds.dims

    def test_open_dataset_with_bounds(self, sample_tropo_product: Path):
        """Should subset by bounds."""
        product = TropoProduct.from_path(sample_tropo_product)
        bounds = (-117.5, 34.25, -117.25, 34.75)
        ds = product.open_dataset(bounds=bounds)

        assert len(ds.latitude) < 50
        assert len(ds.longitude) < 50

    def test_get_total_delay(self, sample_tropo_product: Path):
        """Should compute total delay."""
        product = TropoProduct.from_path(sample_tropo_product)
        total = product.get_total_delay()

        assert total.name == "zenith_total_delay"
        assert "height" in total.dims
        assert "latitude" in total.dims
        assert "longitude" in total.dims


class TestTropoProductMatching:
    """Tests for date matching."""

    def test_matches_date_within_window(self):
        """Should match dates within window."""
        product = TropoProduct(
            path=Path("test.nc"),
            date=datetime(2022, 1, 11, 0, 0),
            production_date=datetime(2022, 1, 11, 12, 0),
            model="HRRR",
            version="1.0",
        )

        target = datetime(2022, 1, 11, 3, 0)

        assert product.matches_date(target, hours=6.0)

    def test_matches_date_outside_window(self):
        """Should not match dates outside window."""
        product = TropoProduct(
            path=Path("test.nc"),
            date=datetime(2022, 1, 11, 0, 0),
            production_date=datetime(2022, 1, 11, 12, 0),
            model="HRRR",
            version="1.0",
        )

        target = datetime(2022, 1, 11, 12, 0)

        assert not product.matches_date(target, hours=6.0)


def test_interpolate_to_dem_surface_is_float32_for_float16_dem():
    """A float16 DEM must not quantise the interpolated delay to ~2 mm steps."""
    import numpy as np
    import rioxarray  # noqa: F401 — registers the .rio accessor
    import xarray as xr

    from cal_disp.product import interpolate_to_dem_surface

    ny, nx = 40, 50
    lat = np.linspace(35.0, 34.0, ny)
    lon = np.linspace(-118.0, -117.0, nx)
    height = np.linspace(0.0, 4000.0, 9)
    # Zenith delay: ~2.4 m at sea level, decreasing with height, plus a
    # gentle horizontal gradient of 0.02 mm per column
    delay = (
        2.4
        - 2.5e-4 * height[:, None, None]
        + 2e-5 * np.arange(nx)[None, None, :]
        + np.zeros((1, ny, 1))
    )
    cube = xr.DataArray(
        delay,
        dims=["height", "latitude", "longitude"],
        coords={"height": height, "latitude": lat, "longitude": lon},
        name="zenith_total_delay",
    ).rio.write_crs("EPSG:4326")

    rng = np.random.default_rng(0)
    dem16 = xr.DataArray(
        rng.uniform(0, 3000, (ny, nx)).astype(np.float16),
        dims=["y", "x"],
        coords={"y": lat, "x": lon},
        name="dem",
        attrs={"units": "m"},
    ).rio.write_crs("EPSG:4326")
    dem16_before = dem16.values.copy()

    out = interpolate_to_dem_surface(cube, dem16)

    assert out.dtype == np.float32
    assert out.rio.crs == dem16.rio.crs
    assert out.shape == dem16.shape
    # The DEM itself is untouched (no write into its float16 buffer)
    np.testing.assert_array_equal(dem16.values, dem16_before)
    # Exact value expected from the (linear) delay model
    expected = 2.4 - 2.5e-4 * dem16.values.astype(np.float64) + 2e-5 * np.arange(nx)
    np.testing.assert_allclose(out.values, expected, atol=2e-6)
    # Not quantised: had the delay been cast into the DEM's float16 buffer
    # (the old behaviour), values near 2 m would sit on ~1-2 mm steps
    quantised = out.values.astype(np.float16).astype(np.float32)
    assert np.abs(quantised - out.values).max() > 4e-4  # >= 0.4 mm error
    assert len(np.unique(out.values)) > 0.9 * out.size
    assert len(np.unique(quantised)) < 0.5 * len(np.unique(out.values))

    # Row-block interpolation (memory bound) gives bit-identical output,
    # including a block size that does not divide the number of rows
    for block_rows in (1, 7, ny + 5):
        blocked = interpolate_to_dem_surface(cube, dem16, block_rows=block_rows)
        np.testing.assert_array_equal(blocked.values, out.values)


def test_repr():
    """Should have readable repr."""
    product = TropoProduct(
        path=Path("test.nc"),
        date=datetime(2022, 1, 11, 0, 0),
        production_date=datetime(2022, 1, 11, 12, 0),
        model="HRRR",
        version="1.0",
    )

    repr_str = repr(product)
    assert "2022-01-11" in repr_str
    assert "HRRR" in repr_str
