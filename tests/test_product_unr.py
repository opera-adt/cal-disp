"""Tests for UnrGrid class."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

from cal_disp.product import UnrGrid


class TestUnrGridBasics:
    """Tests for basic UnrGrid functionality."""

    def test_init(self, sample_unr_data: tuple[Path, Path]):
        """Should initialize with paths."""
        lookup, tenv8_dir = sample_unr_data

        grid = UnrGrid(
            lookup_table=lookup,
            data_dir=tenv8_dir,
            frame_id=8882,
        )

        assert grid.frame_id == 8882
        assert grid.lookup_table.exists()
        assert grid.data_dir.exists()

    def test_missing_lookup(self, tmp_path: Path):
        """Should reject missing lookup table."""
        with pytest.raises(FileNotFoundError, match="Lookup table not found"):
            UnrGrid(
                lookup_table=tmp_path / "nonexistent.txt",
                data_dir=tmp_path,
            )

    def test_missing_data_dir(self, tmp_path: Path):
        """Should reject missing data directory."""
        lookup = tmp_path / "lookup.txt"
        lookup.touch()

        with pytest.raises(FileNotFoundError, match="Data directory not found"):
            UnrGrid(
                lookup_table=lookup,
                data_dir=tmp_path / "nonexistent",
            )


class TestUnrGridDataAccess:
    """Tests for data access methods."""

    def test_to_dataframe(self, sample_unr_data: tuple[Path, Path]):
        """Should construct dataframe."""
        lookup, tenv8_dir = sample_unr_data
        grid = UnrGrid(lookup_table=lookup, data_dir=tenv8_dir)

        df = grid.to_dataframe()

        assert "grid_point" in df.columns
        assert "lon" in df.columns
        assert "lat" in df.columns
        assert "date" in df.columns
        assert "east" in df.columns
        assert "north" in df.columns
        assert "up" in df.columns

    def test_dataframe_units(self, sample_unr_data: tuple[Path, Path]):
        """Should convert to meters."""
        lookup, tenv8_dir = sample_unr_data
        grid = UnrGrid(lookup_table=lookup, data_dir=tenv8_dir)

        df = grid.to_dataframe()

        # Original tenv8 has 1.2mm, should be 0.0012m
        assert df["east"].max() == pytest.approx(0.0024, abs=0.0001)

    def test_get_grid_points(self, sample_unr_data: tuple[Path, Path]):
        """Should get unique grid points."""
        lookup, tenv8_dir = sample_unr_data
        grid = UnrGrid(lookup_table=lookup, data_dir=tenv8_dir)

        points = grid.get_grid_points()

        assert len(points) == 4
        assert all(col in points.columns for col in ["grid_point", "lon", "lat"])

    def test_get_bounds(self, sample_unr_data: tuple[Path, Path]):
        """Should get spatial bounds."""
        lookup, tenv8_dir = sample_unr_data
        grid = UnrGrid(lookup_table=lookup, data_dir=tenv8_dir)

        west, south, east, north = grid.get_bounds()

        assert west == pytest.approx(-118.0)
        assert east == pytest.approx(-117.9)
        assert south == pytest.approx(34.0)
        assert north == pytest.approx(34.1)

    def test_get_time_range(self, sample_unr_data: tuple[Path, Path]):
        """Should get temporal range."""
        lookup, tenv8_dir = sample_unr_data
        grid = UnrGrid(lookup_table=lookup, data_dir=tenv8_dir)

        start, end = grid.get_time_range()

        assert isinstance(start, datetime)
        assert isinstance(end, datetime)
        assert start < end


def test_repr(sample_unr_data: tuple[Path, Path]):
    """Should have readable repr."""
    lookup, tenv8_dir = sample_unr_data
    grid = UnrGrid(lookup_table=lookup, data_dir=tenv8_dir, frame_id=8882)

    repr_str = repr(grid)

    # Just check basic info is present
    assert "8882" in repr_str
    assert "UnrGrid" in repr_str
