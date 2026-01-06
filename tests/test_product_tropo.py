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
