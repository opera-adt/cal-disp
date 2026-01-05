from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import pandas as pd

from cal_disp.download._stage_unr import load_lookup_table
from cal_disp.product._utils import decimal_year_to_datetime


@dataclass
class UnrGrid:
    """UNR GNSS grid data from lookup table and tenv8 files.

    Represents gridded GNSS displacement timeseries from University of Nevada Reno.
    Combines geographic coordinates from lookup table with displacement data
    from .tenv8 files.

    Parameters
    ----------
    lookup_table : Path
        Path to grid_latlon_lookup.txt file.
    data_dir : Path
        Directory containing .tenv8 files.
    frame_id : int or None, optional
        OPERA frame identifier. Default is None.

    Examples
    --------
    >>> grid = UnrGrid(
    ...     lookup_table=Path("data/grid_latlon_lookup_v0.2.txt"),
    ...     data_dir=Path("data/tenv8_files"),
    ...     frame_id=8882
    ... )
    >>> df = grid.to_dataframe()
    >>> df.columns
    ['grid_point', 'lon', 'lat', 'date', 'east', 'north', 'up', ...]

    """

    lookup_table: Path
    data_dir: Path
    frame_id: int | None = None

    _TENV8_PATTERN = re.compile(r"(\d{6})_[A-Z0-9]+\.tenv8")

    def __post_init__(self) -> None:
        """Validate paths after construction."""
        self.lookup_table = Path(self.lookup_table)
        self.data_dir = Path(self.data_dir)

        if not self.lookup_table.exists():
            raise FileNotFoundError(f"Lookup table not found: {self.lookup_table}")
        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {self.data_dir}")

        self._df_cache = None

    def _load_tenv8(self, path: Path) -> pd.DataFrame:
        """Load a single tenv8 file.

        Parameters
        ----------
        path : Path
            Path to .tenv8 file.

        Returns
        -------
        pd.DataFrame
            Timeseries data with grid_point column added.

        """
        match = self._TENV8_PATTERN.search(path.name)
        if not match:
            raise ValueError(f"Cannot parse grid_point from filename: {path.name}")

        grid_point = int(match.group(1))

        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            names=[
                "decimal_year",
                "east",
                "north",
                "up",
                "sigma_east",
                "sigma_north",
                "sigma_up",
                "rapid_flag",
            ],
        )
        df["grid_point"] = grid_point
        return df

    def to_dataframe(self, use_cache: bool = True) -> pd.DataFrame:
        """Construct grid dataframe from lookup table and tenv8 files.

        Combines geographic coordinates with displacement timeseries data.
        Converts decimal years to datetime and displacement values from
        millimeters to meters.

        Parameters
        ----------
        use_cache : bool, optional
            If True, cache result for subsequent calls. Default is True.

        Returns
        -------
        pd.DataFrame
            Combined dataframe with columns:
            - grid_point: Grid point ID
            - lon, lat: Geographic coordinates (normalized to [-180, 180])
            - date: Observation datetime
            - east, north, up: Displacement components (meters)
            - sigma_east, sigma_north, sigma_up: Uncertainties (meters)
            - corr_en, corr_eu, corr_nu: Correlation coefficients (placeholder)
            - rapid_flag: Rapid solution flag

        Raises
        ------
        FileNotFoundError
            If no tenv8 files found in data directory.

        Notes
        -----
        UNR data is in millimeters; this method converts to meters.
        Correlation coefficients are set to 0.0 (not provided in .tenv8 format).

        Examples
        --------
        >>> grid = UnrGrid(lookup_table=Path("lookup.txt"), data_dir=Path("data"))
        >>> df = grid.to_dataframe()
        >>> df['east'].max()  # in meters
        0.045

        """
        if use_cache and self._df_cache is not None:
            return self._df_cache

        lookup = load_lookup_table(self.lookup_table, normalize_longitude=True)

        tenv8_files = sorted(self.data_dir.glob("*.tenv8"))
        if not tenv8_files:
            raise FileNotFoundError(f"No .tenv8 files found in {self.data_dir}")

        timeseries = pd.concat(
            [self._load_tenv8(path) for path in tenv8_files],
            ignore_index=True,
        )

        timeseries["date"] = timeseries["decimal_year"].apply(decimal_year_to_datetime)

        # Placeholder correlation values (not in .tenv8 format)
        timeseries["corr_en"] = 0.0
        timeseries["corr_eu"] = 0.0
        timeseries["corr_nu"] = 0.0

        result = timeseries.merge(lookup.reset_index(), on="grid_point", how="left")

        # Convert from millimeters to meters
        for col in ["east", "north", "up", "sigma_east", "sigma_north", "sigma_up"]:
            result[col] /= 1000.0

        result = result[
            [
                "grid_point",
                "lon",
                "lat",
                "date",
                "east",
                "north",
                "up",
                "sigma_east",
                "sigma_north",
                "sigma_up",
                "corr_en",
                "corr_eu",
                "corr_nu",
                "rapid_flag",
            ]
        ]

        result.attrs["units"] = "meters"

        if use_cache:
            self._df_cache = result

        return result

    def to_geodataframe(self) -> gpd.GeoDataFrame:
        """Convert to GeoDataFrame with point geometries.

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame with point geometries at each grid location.

        Examples
        --------
        >>> grid = UnrGrid(lookup_table=Path("lookup.txt"), data_dir=Path("data"))
        >>> gdf = grid.to_geodataframe()
        >>> gdf.crs
        'EPSG:4326'

        """
        df = self.to_dataframe()
        return gpd.GeoDataFrame(
            df,
            geometry=gpd.points_from_xy(x=df.lon, y=df.lat),
            crs="EPSG:4326",
        )

    def get_grid_points(self) -> pd.DataFrame:
        """Get unique grid points with their coordinates.

        Returns
        -------
        pd.DataFrame
            DataFrame with one row per grid point: grid_point, lon, lat.

        Examples
        --------
        >>> grid = UnrGrid(lookup_table=Path("lookup.txt"), data_dir=Path("data"))
        >>> points = grid.get_grid_points()
        >>> len(points)
        450

        """
        df = self.to_dataframe()
        return df[["grid_point", "lon", "lat"]].drop_duplicates().reset_index(drop=True)

    def get_bounds(self) -> tuple[float, float, float, float]:
        """Get spatial bounds of grid.

        Returns
        -------
        tuple[float, float, float, float]
            (west, south, east, north) bounds.

        Examples
        --------
        >>> grid.get_bounds()
        (-120.5, 32.1, -115.2, 35.8)

        """
        points = self.get_grid_points()
        return (
            points.lon.min(),
            points.lat.min(),
            points.lon.max(),
            points.lat.max(),
        )

    def get_time_range(self) -> tuple[datetime, datetime]:
        """Get temporal range of observations.

        Returns
        -------
        tuple[datetime, datetime]
            (min_date, max_date).

        Examples
        --------
        >>> grid.get_time_range()
        (datetime(2014, 7, 1), datetime(2024, 10, 15))

        """
        df = self.to_dataframe()
        return (df.date.min(), df.date.max())

    def __repr__(self) -> str:
        """Return a string representation."""
        frame_str = f"frame={self.frame_id}" if self.frame_id else "no frame"
        try:
            n_points = len(self.get_grid_points())
            start, end = self.get_time_range()
            return (
                f"UnrGrid({frame_str}, points={n_points}, "
                f"dates={start.date()}-{end.date()})"
            )
        except Exception:
            return f"UnrGrid({frame_str})"
