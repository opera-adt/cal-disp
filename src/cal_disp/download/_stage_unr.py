from __future__ import annotations

import logging
from functools import partial
from pathlib import Path
from typing import Literal, Sequence

import geopandas as gpd
import pandas as pd
import requests
from opera_utils import get_frame_geojson
from requests.adapters import HTTPAdapter
from shapely.geometry import box
from tqdm.contrib.concurrent import thread_map
from urllib3.util.retry import Retry

__all__ = [
    "create_session",
    "download_lookup_table",
    "load_lookup_table",
    "download_grid_file",
    "download_grid_files",
    "download_unr_grid",
]

logger = logging.getLogger(__name__)

# Constants
VALID_VERSIONS = {"0.1", "0.2"}
DEFAULT_VERSION: Literal["0.1", "0.2"] = "0.2"

LOOKUP_URL = (
    "https://geodesy.unr.edu/grid_timeseries/Version{version}/grid_latlon_lookup.txt"
)
GRID_BASE_URL = (
    "https://geodesy.unr.edu/grid_timeseries/Version{version}/time_variable_gridded"
)

# Type aliases
PlateType = Literal["NA", "PA", "IGS14", "IGS20"]
VersionType = Literal["0.1", "0.2"]


def create_session(retries: int = 5, backoff: float = 1.0) -> requests.Session:
    """Create a requests session with retry logic.

    Parameters
    ----------
    retries : int, optional
        Number of retry attempts. Default is 5.
    backoff : float, optional
        Backoff factor between retries. Default is 1.0.

    Returns
    -------
    requests.Session
        Configured session with retry adapter.

    """
    session = requests.Session()
    retry_strategy = Retry(
        total=retries,
        backoff_factor=backoff,
        status_forcelist=[502, 503, 504],
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("https://", adapter)
    return session


def download_lookup_table(
    output_dir: Path,
    version: VersionType = DEFAULT_VERSION,
    session: requests.Session | None = None,
) -> Path:
    """Download the UNR grid latitude/longitude lookup table.

    This table maps grid point IDs to geographic coordinates. The file
    is saved in the original space-separated format.

    Parameters
    ----------
    output_dir : Path
        Directory where lookup file will be saved.
    version : {"0.1", "0.2"}, optional
        UNR data version. Default is "0.2".
    session : requests.Session or None, optional
        Session with retry logic. If None, a new session is created.

    Returns
    -------
    Path
        Path to the downloaded lookup file.

    Raises
    ------
    ValueError
        If version is not supported.
    requests.HTTPError
        If download fails.

    Notes
    -----
    The file format is space-separated with three columns:
    grid_point longitude latitude (no header).

    Version 0.2 uses longitude range [0, 360] in the original file.
    This function preserves the original format without modification.

    Examples
    --------
    >>> from pathlib import Path
    >>> lookup_path = download_lookup_table(Path("data"))
    >>> df = load_lookup_table(lookup_path)

    """
    if version not in VALID_VERSIONS:
        msg = f"Version must be one of {VALID_VERSIONS}, got '{version}'"
        raise ValueError(msg)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"grid_latlon_lookup_v{version}.txt"

    # Skip if already downloaded
    if output_path.exists():
        logger.debug(f"Lookup table already exists: {output_path}")
        return output_path

    url = LOOKUP_URL.format(version=version)
    logger.info(f"Downloading lookup table from {url}")

    if session is None:
        session = create_session()

    response = session.get(url)
    response.raise_for_status()

    output_path.write_bytes(response.content)
    logger.info(f"Saved lookup table to {output_path}")

    return output_path


def load_lookup_table(path: Path, normalize_longitude: bool = True) -> pd.DataFrame:
    """Load a lookup table file into a DataFrame.

    Parameters
    ----------
    path : Path
        Path to lookup table file.
    normalize_longitude : bool, optional
        If True, convert longitude from [0, 360] to [-180, 180].
        Default is True.

    Returns
    -------
    pd.DataFrame
        Lookup table with columns: grid_point (index), lon, lat, alt.

    Examples
    --------
    >>> lookup_path = download_lookup_table(Path("data"))
    >>> df = load_lookup_table(lookup_path)
    >>> lat, lon = df.loc[123456, ['lat', 'lon']]

    """
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=["grid_point", "lon", "lat"],
    )

    # Normalize longitude to [-180, 180] if requested
    if normalize_longitude:
        df["lon"] = ((df["lon"] + 180) % 360) - 180

    # Grid points don't have altitude info
    df["alt"] = 0.0

    return df.set_index("grid_point")


def download_grid_file(
    grid_id: int,
    output_dir: Path,
    plate: PlateType = "IGS20",
    version: VersionType = DEFAULT_VERSION,
    session: requests.Session | None = None,
) -> Path:
    r"""Download a single grid point timeseries file.

    Downloads a .tenv8 file containing displacement timeseries data
    for the specified grid point.

    Parameters
    ----------
    grid_id : int
        Grid point identifier (e.g., 123456).
    output_dir : Path
        Directory where file will be saved.
    plate : {"NA", "PA", "IGS14", "IGS20"}, optional
        Reference plate for the data. Default is "IGS20".
    version : {"0.1", "0.2"}, optional
        UNR data version. Default is "0.2".
    session : requests.Session or None, optional
        Session with retry logic. If None, a new session is created.

    Returns
    -------
    Path
        Path to the downloaded file.

    Raises
    ------
    ValueError
        If version is not supported.
    requests.HTTPError
        If download fails.

    Notes
    -----
    IGS14 plate is not available in version 0.2. The function automatically
    uses IGS20 instead if IGS14 is requested with version 0.2.

    The .tenv8 format contains columns: decimal_year, east, north, up,
    sigma_east, sigma_north, sigma_up, rapid_flag.

    Examples
    --------
    >>> from pathlib import Path
    >>> output = download_grid_file(123456, Path("data"))
    >>> df = pd.read_csv(output, sep=r"\s+", header=None)

    """
    if version not in VALID_VERSIONS:
        msg = f"Version must be one of {VALID_VERSIONS}, got '{version}'"
        raise ValueError(msg)

    # Handle IGS14/IGS20 plate compatibility
    if plate == "IGS14" and version == "0.2":
        plate = "IGS20"

    # Build URL and output path
    filename = f"{plate}/{grid_id:06d}_{plate}.tenv8"
    url = f"{GRID_BASE_URL.format(version=version)}/{filename}"

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{grid_id:06d}_{plate}.tenv8"

    # Skip if already downloaded
    if output_path.exists():
        return output_path

    # Download
    if session is None:
        session = create_session()

    response = session.get(url)
    response.raise_for_status()

    output_path.write_bytes(response.content)

    return output_path


def download_grid_files(
    grid_ids: Sequence[int],
    output_dir: Path,
    plate: PlateType = "IGS20",
    version: VersionType = DEFAULT_VERSION,
    max_workers: int = 4,
) -> list[Path]:
    """Download multiple grid files in parallel.

    Parameters
    ----------
    grid_ids : Sequence[int]
        List of grid point IDs to download.
    output_dir : Path
        Directory where files will be saved.
    plate : {"NA", "PA", "IGS14", "IGS20"}, optional
        Reference plate. Default is "IGS20".
    version : {"0.1", "0.2"}, optional
        UNR data version. Default is "0.2".
    max_workers : int, optional
        Number of parallel download threads. Default is 4.

    Returns
    -------
    list[Path]
        Paths to downloaded files.

    Examples
    --------
    >>> grid_ids = [123456, 123457, 123458]
    >>> paths = download_grid_files(grid_ids, Path("data"), max_workers=8)

    """
    logger.info(f"Downloading {len(grid_ids)} grid files with {max_workers} workers")

    # Create shared session for all downloads
    session = create_session()

    # Fix constant arguments
    worker = partial(
        download_grid_file,
        output_dir=output_dir,
        plate=plate,
        version=version,
        session=session,
    )

    # Download in parallel with progress bar
    paths = thread_map(
        worker,
        grid_ids,
        max_workers=max_workers,
        desc="Downloading grid files",
    )

    return list(paths)


def download_unr_grid(
    frame_id: int,
    output_dir: Path,
    margin_deg: float = 0.5,
    plate: PlateType = "IGS20",
    version: VersionType = DEFAULT_VERSION,
    max_workers: int = 4,
) -> Path:
    """Download UNR gridded GNSS timeseries for a given frame.

    Downloads .tenv8 files for all grid points within the frame bounds.
    Use UnrGrid to load the downloaded data.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    output_dir : Path
        Output directory for downloaded data.
    margin_deg : float, optional
        Margin in degrees to expand frame bounding box. Default is 0.5.
    plate : {"NA", "PA", "IGS14", "IGS20"}, optional
        Reference plate. Default is "IGS20".
    version : {"0.1", "0.2"}, optional
        UNR grid version. Default is "0.2".
    max_workers : int, optional
        Number of parallel download threads. Default is 4.

    Returns
    -------
    Path
        Directory containing downloaded .tenv8 files and lookup table.

    Examples
    --------
    >>> data_dir = download_unr_grid(8882, Path("data"))
    >>> # Load with UnrGrid
    >>> from .grid import UnrGrid
    >>> grid = UnrGrid(
    ...     lookup_table=data_dir / "grid_latlon_lookup_v0.2.txt",
    ...     data_dir=data_dir,
    ...     frame_id=8882
    ... )
    >>> df = grid.to_dataframe()

    """
    logger.info(
        f"Downloading UNR grid for frame {frame_id} "
        f"(plate={plate}, version={version}, margin={margin_deg}°)"
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create session with retry logic
    session = create_session()

    # Download lookup table
    lookup_file = download_lookup_table(
        output_dir=output_dir,
        version=version,
        session=session,
    )

    # Load lookup and convert to GeoDataFrame
    if version == "0.2":
        normalize_longitude = True
    else:
        normalize_longitude = False

    lookup = load_lookup_table(lookup_file, normalize_longitude=normalize_longitude)
    grid_gdf = gpd.GeoDataFrame(
        lookup,
        geometry=gpd.points_from_xy(x=lookup.lon, y=lookup.lat),
        crs="EPSG:4326",
    )

    # Get frame bounds and expand by margin
    frame_gdf = get_frame_geojson([frame_id], as_geodataframe=True)
    west, south, east, north = frame_gdf.bounds.values[0]
    bounds_poly = box(
        west - margin_deg,
        south - margin_deg,
        east + margin_deg,
        north + margin_deg,
    )

    # Filter grid points to expanded frame bounds
    grid_gdf = grid_gdf.clip(bounds_poly)
    grid_ids = grid_gdf.index.tolist()

    logger.info(f"Found {len(grid_ids)} grid points within frame bounds")

    # Download .tenv8 files for all grid points in parallel
    download_grid_files(
        grid_ids=grid_ids,
        output_dir=output_dir,
        plate=plate,
        version=version,
        max_workers=max_workers,
    )

    logger.info(f"Download complete. Data saved to {output_dir}")

    return output_dir
