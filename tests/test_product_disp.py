"""Tests for DispProduct class."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

from cal_disp.product import DispProduct


class TestDispProductParsing:
    """Tests for filename parsing."""

    def test_from_path_s1(self):
        """Should parse S1 filename correctly."""
        filename = (
            "OPERA_L3_DISP-S1_IW_F08882_VV_"
            "20220111T002651Z_20220722T002657Z_"
            "v1.0_20230101T000000Z.nc"
        )
        product = DispProduct.from_path(filename)

        assert product.frame_id == 8882
        assert product.sensor == "S1"
        assert product.mode == "IW"
        assert product.polarization == "VV"
        assert product.version == "1.0"
        assert product.primary_date == datetime(2022, 1, 11, 0, 26, 51)
        assert product.secondary_date == datetime(2022, 7, 22, 0, 26, 57)

    def test_from_path_invalid_filename(self):
        """Should reject invalid filename."""
        with pytest.raises(ValueError, match="does not match"):
            DispProduct.from_path("invalid_filename.nc")

    def test_validates_frame_id(self):
        """Should reject invalid frame ID."""
        with pytest.raises(ValueError, match="must be positive"):
            DispProduct(
                path=Path("test.nc"),
                frame_id=0,
                primary_date=datetime(2022, 1, 1),
                secondary_date=datetime(2022, 2, 1),
                sensor="S1",
                polarization="VV",
                version="1.0",
                production_date=datetime(2023, 1, 1),
            )

    def test_validates_date_order(self):
        """Should reject secondary before primary."""
        with pytest.raises(ValueError, match="must be after"):
            DispProduct(
                path=Path("test.nc"),
                frame_id=8882,
                primary_date=datetime(2022, 2, 1),
                secondary_date=datetime(2022, 1, 1),
                sensor="S1",
                polarization="VV",
                version="1.0",
                production_date=datetime(2023, 1, 1),
            )

    def test_validates_polarization(self):
        """Should reject invalid polarization."""
        with pytest.raises(ValueError, match="Invalid polarization"):
            DispProduct(
                path=Path("test.nc"),
                frame_id=8882,
                primary_date=datetime(2022, 1, 1),
                secondary_date=datetime(2022, 2, 1),
                sensor="S1",
                polarization="XX",
                version="1.0",
                production_date=datetime(2023, 1, 1),
            )


class TestDispProductProperties:
    """Tests for product properties."""

    def test_baseline_days(self):
        """Should calculate temporal baseline."""
        product = DispProduct(
            path=Path("test.nc"),
            frame_id=8882,
            primary_date=datetime(2022, 1, 1),
            secondary_date=datetime(2022, 1, 13),
            sensor="S1",
            polarization="VV",
            version="1.0",
            production_date=datetime(2023, 1, 1),
        )

        assert product.baseline_days == 12

    def test_filename_property(self):
        """Should return filename from path."""
        product = DispProduct(
            path=Path("/data/test.nc"),
            frame_id=8882,
            primary_date=datetime(2022, 1, 1),
            secondary_date=datetime(2022, 2, 1),
            sensor="S1",
            polarization="VV",
            version="1.0",
            production_date=datetime(2023, 1, 1),
        )

        assert product.filename == "test.nc"

    def test_exists_false(self):
        """Should return False for non-existent file."""
        product = DispProduct(
            path=Path("/nonexistent/test.nc"),
            frame_id=8882,
            primary_date=datetime(2022, 1, 1),
            secondary_date=datetime(2022, 2, 1),
            sensor="S1",
            polarization="VV",
            version="1.0",
            production_date=datetime(2023, 1, 1),
        )

        assert not product.exists


class TestDispProductDataAccess:
    """Tests for data access methods."""

    def test_open_dataset(self, sample_disp_product_with_corrections: Path):
        """Should open main dataset."""
        product = DispProduct.from_path(sample_disp_product_with_corrections)
        ds = product.open_dataset()

        assert "displacement" in ds
        assert "temporal_coherence" in ds
        assert "spatial_ref" in ds

    def test_open_corrections(self, sample_disp_product_with_corrections: Path):
        """Should open corrections group."""
        product = DispProduct.from_path(sample_disp_product_with_corrections)
        ds = product.open_corrections()

        assert "ionospheric_delay" in ds
        assert "solid_earth_tide" in ds
        assert "reference_point" in ds

    def test_get_reference_point_index(
        self, sample_disp_product_with_corrections: Path
    ):
        """Should get reference point indices."""
        product = DispProduct.from_path(sample_disp_product_with_corrections)
        row, col = product.get_reference_point_index()

        assert row == 25
        assert col == 25

    def test_get_reference_point_latlon(
        self, sample_disp_product_with_corrections: Path
    ):
        """Should get reference point coordinates."""
        product = DispProduct.from_path(sample_disp_product_with_corrections)
        lat, lon = product.get_reference_point_latlon()

        assert lat == pytest.approx(34.5)
        assert lon == pytest.approx(-118.0)

    def test_get_epsg(self, sample_disp_product_with_corrections: Path):
        """Should extract EPSG code."""
        product = DispProduct.from_path(sample_disp_product_with_corrections)
        epsg = product.get_epsg()

        assert epsg == 32615  # UTM zone 15N

    def test_get_bounds(self, sample_disp_product_with_corrections: Path):
        """Should get native bounds."""
        product = DispProduct.from_path(sample_disp_product_with_corrections)
        bounds = product.get_bounds()

        assert "left" in bounds
        assert "bottom" in bounds
        assert "right" in bounds
        assert "top" in bounds
        assert bounds["left"] < bounds["right"]
        assert bounds["bottom"] < bounds["top"]


def test_repr():
    """Should have readable repr."""
    product = DispProduct(
        path=Path("test.nc"),
        frame_id=8882,
        primary_date=datetime(2022, 1, 1),
        secondary_date=datetime(2022, 1, 13),
        sensor="S1",
        polarization="VV",
        version="1.0",
        production_date=datetime(2023, 1, 1),
    )

    repr_str = repr(product)
    assert "8882" in repr_str
    assert "VV" in repr_str
