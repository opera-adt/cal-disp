import re
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
import pandas as pd
import pyarrow.parquet as pq


@dataclass
class UnrGrid:
    """UNR GNSS grid data.

    Represents gridded GNSS velocity data from University of Nevada Reno.
    Data is stored as parquet with point geometries and metadata.

    Examples
    --------
    >>> # Load from path (frame_id parsed if in filename)
    >>> grid = UnrGrid.from_path("unr_grid_frame8882.parquet")
    >>> grid.frame_id
    8882

    >>> # Load GeoDataFrame
    >>> gdf = grid.load()
    >>> gdf.columns
    ['lon', 'lat', 'east', 'north', 'up', 'geometry', ...]

    >>> # Get metadata
    >>> meta = grid.get_metadata()
    >>> meta['source']
    'UNR'

    """

    path: Path
    frame_id: int | None = None

    # Optional pattern to extract frame_id from filename
    _PATTERN = re.compile(r"frame[\s_-]?(\d+)", re.IGNORECASE)

    def __post_init__(self) -> None:
        """Validate grid after construction."""
        self.path = Path(self.path)

    @classmethod
    def from_path(cls, path: Path | str, frame_id: int | None = None) -> "UnrGrid":
        """Create UnrGrid from parquet file path.

        Parameters
        ----------
        path : Path or str
            Path to UNR parquet file.
        frame_id : int or None, optional
            Frame ID. If None, attempts to parse from filename.
            Default is None.

        Returns
        -------
        UnrGrid
            Grid instance.

        Examples
        --------
        >>> # Frame ID from filename
        >>> grid = UnrGrid.from_path("unr_grid_frame8882.parquet")
        >>> grid.frame_id
        8882

        >>> # Explicit frame ID
        >>> grid = UnrGrid.from_path("custom_unr_data.parquet", frame_id=8882)
        >>> grid.frame_id
        8882

        >>> # No frame ID
        >>> grid = UnrGrid.from_path("unr_data.parquet")
        >>> grid.frame_id is None
        True

        """
        path = Path(path)

        # Try to parse frame_id from filename if not provided
        if frame_id is None:
            match = cls._PATTERN.search(path.name)
            if match:
                frame_id = int(match.group(1))

        return cls(path=path, frame_id=frame_id)

    def load(self) -> gpd.GeoDataFrame:
        """Load UNR grid as GeoDataFrame.

        Returns
        -------
        gpd.GeoDataFrame
            GeoDataFrame with point geometries and velocity data.

        Raises
        ------
        FileNotFoundError
            If parquet file does not exist.

        Examples
        --------
        >>> grid = UnrGrid.from_path("unr_grid_frame8882.parquet")
        >>> gdf = grid.load()
        >>> gdf.crs
        'EPSG:4326'
        >>> gdf[['lon', 'lat', 'east', 'north', 'up']].head()

        """
        if not self.path.exists():
            raise FileNotFoundError(f"UNR grid file not found: {self.path}")

        # Load parquet as DataFrame
        df = pd.read_parquet(self.path)

        # Create GeoDataFrame with point geometries
        gdf = gpd.GeoDataFrame(
            df,
            geometry=gpd.points_from_xy(x=df.lon, y=df.lat),
            crs="EPSG:4326",
        )

        return gdf

    def get_metadata(self) -> dict[str, str]:
        """Extract metadata from parquet file.

        Returns
        -------
        dict[str, str]
            Metadata dictionary.

        Examples
        --------
        >>> grid = UnrGrid.from_path("unr_grid_frame8882.parquet")
        >>> meta = grid.get_metadata()
        >>> meta.keys()
        dict_keys(['source', 'date_created', 'frame_id', ...])

        """
        if not self.path.exists():
            raise FileNotFoundError(f"UNR grid file not found: {self.path}")

        meta = pq.read_metadata(self.path).metadata

        if meta is None:
            return {}

        metadata_dict = {k.decode(): v.decode() for k, v in meta.items()}

        return metadata_dict

    def to_dataframe(self) -> pd.DataFrame:
        """Load as regular DataFrame without geometry.

        Returns
        -------
        pd.DataFrame
            DataFrame with lon, lat, and velocity columns.

        """
        if not self.path.exists():
            raise FileNotFoundError(f"UNR grid file not found: {self.path}")

        return pd.read_parquet(self.path)

    def get_bounds(self) -> dict[str, float]:
        """Get spatial bounds of grid.

        Returns
        -------
        dict[str, float]
            Dictionary with keys: west, south, east, north.

        """
        gdf = self.load()
        bounds = gdf.total_bounds  # (minx, miny, maxx, maxy)

        return {
            "west": bounds[0],
            "south": bounds[1],
            "east": bounds[2],
            "north": bounds[3],
        }

    def get_grid_count(self) -> int:
        """Get number of GNSS points in grid.

        Returns
        -------
        int
            Number of stations.

        """
        df = self.to_dataframe()
        grid_points = df.groupby("id", as_index=False).first()
        return len(grid_points)

    @property
    def filename(self) -> str:
        """Grid filename."""
        return self.path.name

    @property
    def exists(self) -> bool:
        """Check if grid file exists."""
        return self.path.exists()

    def __repr__(self) -> str:
        """Return a string representation."""
        frame_str = f"frame={self.frame_id}" if self.frame_id else "frame=None"
        return (
            f"UnrGrid({frame_str},"
            f" points={self.get_grid_count() if self.exists else '?'})"
        )
