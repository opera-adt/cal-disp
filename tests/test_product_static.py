"""Tests for StaticLayer class."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from cal_disp.product import StaticLayer


class TestStaticLayerParsing:
    """Tests for filename parsing."""

    def test_from_path_dem(self):
        """Should parse DEM filename."""
        filename = "OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_dem.tif"
        layer = StaticLayer.from_path(filename)

        assert layer.frame_id == 8882
        assert layer.satellite == "S1A"
        assert layer.version == "1.0"
        assert layer.layer_type == "dem"
        assert layer.reference_date == datetime(2014, 4, 3)

    def test_from_path_los(self):
        """Should parse LOS filename."""
        filename = (
            "OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_line_of_sight_enu.tif"
        )
        layer = StaticLayer.from_path(filename)

        assert layer.layer_type == "line_of_sight_enu"

    def test_from_path_invalid(self):
        """Should reject invalid filename."""
        with pytest.raises(ValueError, match="does not match"):
            StaticLayer.from_path("invalid_file.tif")


class TestStaticLayerDataAccess:
    """Tests for data access methods."""

    def test_read_single_band(self, sample_static_dem: Path):
        """Should read DEM data."""
        layer = StaticLayer.from_path(sample_static_dem)
        data = layer.read(band=1)

        assert data.shape == (100, 100)
        assert data.dtype == np.float32

    def test_read_multiple_bands(self, sample_static_los: Path):
        """Should read all LOS bands."""
        layer = StaticLayer.from_path(sample_static_los)
        bands = layer.read_bands()

        assert len(bands) == 3
        assert all(b.shape == (100, 100) for b in bands)

    def test_num_bands_property(self, sample_static_los: Path):
        """Should report correct band count."""
        layer = StaticLayer.from_path(sample_static_los)

        assert layer.num_bands == 3

    def test_to_dataset_dem(self, sample_static_dem: Path):
        """Should convert DEM to xarray."""
        layer = StaticLayer.from_path(sample_static_dem)
        ds = layer.to_dataset()

        assert "dem" in ds
        assert "x" in ds.coords
        assert "y" in ds.coords
        assert ds["dem"].shape == (100, 100)

    def test_to_dataset_los(self, sample_static_los: Path):
        """Should convert LOS to xarray with components."""
        layer = StaticLayer.from_path(sample_static_los)
        ds = layer.to_dataset()

        assert "los_east" in ds
        assert "los_north" in ds
        assert "los_up" in ds

    def test_compute_incidence_angle(self, sample_static_los: Path):
        """Should compute incidence angle from LOS."""
        layer = StaticLayer.from_path(sample_static_los)
        inc_angle = layer.compute_incidence_angle()

        assert inc_angle.shape == (100, 100)
        assert np.all((inc_angle >= 0) & (inc_angle <= 90))

    def test_compute_incidence_angle_wrong_layer(self, sample_static_dem: Path):
        """Should reject non-LOS layers."""
        layer = StaticLayer.from_path(sample_static_dem)

        with pytest.raises(ValueError, match="only valid for line_of_sight"):
            layer.compute_incidence_angle()


class TestStaticLayerMetadata:
    """Tests for metadata access."""

    def test_get_shape(self, sample_static_dem: Path):
        """Should get array shape."""
        layer = StaticLayer.from_path(sample_static_dem)
        height, width = layer.get_shape()

        assert height == 100
        assert width == 100

    def test_get_bounds(self, sample_static_dem: Path):
        """Should get native bounds."""
        layer = StaticLayer.from_path(sample_static_dem)
        bounds = layer.get_bounds()

        assert "left" in bounds
        assert "bottom" in bounds
        assert "right" in bounds
        assert "top" in bounds

    def test_get_epsg(self, sample_static_dem: Path):
        """Should get EPSG code."""
        layer = StaticLayer.from_path(sample_static_dem)
        epsg = layer.get_epsg()

        assert epsg == 4326

    def test_filename_property(self, sample_static_dem: Path):
        """Should return filename."""
        layer = StaticLayer.from_path(sample_static_dem)

        assert layer.filename == sample_static_dem.name

    def test_exists_property(self, sample_static_dem: Path):
        """Should check file existence."""
        layer = StaticLayer.from_path(sample_static_dem)

        assert layer.exists


def test_repr(sample_static_dem: Path):
    """Should have readable repr."""
    layer = StaticLayer.from_path(sample_static_dem)
    repr_str = repr(layer)

    assert "8882" in repr_str
    assert "dem" in repr_str
