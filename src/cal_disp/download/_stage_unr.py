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
from tqdm.contrib.concurrent import thread_map
from urllib3.util.retry import Retry

from ._errors import DownloadError

__all__ = [
    "DEFAULT_TIMEOUT",
    "GRID_BASE_URLS",
    "DownloadError",
    "create_session",
    "download_file",
    "download_lookup_table",
    "load_lookup_table",
    "download_grid_file",
    "download_grid_files",
    "download_unr_grid",
    "grid_file_name",
]

logger = logging.getLogger(__name__)

# (connect, read) timeout in seconds for every request. Without a timeout a
# stalled server keeps a PGE stage waiting forever (requests has no default).
DEFAULT_TIMEOUT: tuple[float, float] = (10.0, 120.0)
# Statuses worth a retry with backoff: rate limiting, transient server errors
RETRY_STATUSES = (429, 500, 502, 503, 504)
# Suffix of the file a download is streamed to before the atomic rename
PART_SUFFIX = ".part"
CHUNK_SIZE = 1 << 16

# Constants
VALID_VERSIONS = {"0.1", "0.2", "0.3"}
DEFAULT_VERSION: Literal["0.1", "0.2", "0.3"] = "0.3"

LOOKUP_URL = (
    "https://geodesy.unr.edu/grid_timeseries/Version{version}/grid_latlon_lookup.txt"
)
# "contsant" (sic) is UNR's actual path spelling. The constant product holds
# precomputed linear rates and exists only for Version0.3 (IGS20 and NA).
GRID_BASE_URLS = {
    "constant": (
        "https://geodesy.unr.edu/grid_timeseries/Version{version}/time_contsant_gridded"
    ),
    "variable": (
        "https://geodesy.unr.edu/grid_timeseries/Version{version}/time_variable_gridded"
    ),
}
CONSTANT_VERSIONS = {"0.3"}
CONSTANT_PLATES = {"IGS20", "NA"}

# Type aliases
PlateType = Literal["NA", "PA", "IGS14", "IGS20"]
VersionType = Literal["0.1", "0.2", "0.3"]
GridType = Literal["constant", "variable"]
DEFAULT_GRID_TYPE: GridType = "constant"


def grid_file_name(grid_id: int, plate: str, grid_type: str) -> str:
    """Local file name of a staged grid point, ``<id>_<plate>_<grid_type>.tenv8``.

    Both UNR products share one remote file name, so the grid type is part of
    the local name. This matches the name Venti's ``download_station`` looks
    for, so Venti reuses staged files instead of downloading them again.
    """
    return f"{grid_id:06d}_{plate}_{grid_type}.tenv8"


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
        status_forcelist=list(RETRY_STATUSES),
        allowed_methods=["GET", "HEAD"],
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def _looks_like_html(head: bytes) -> bool:
    """Whether the first bytes of a body are an HTML document."""
    start = head.lstrip()[:32].lower()
    return start.startswith((b"<!doctype html", b"<html", b"<head", b"<body"))


def download_file(
    url: str,
    output_path: Path,
    session: requests.Session | None = None,
    timeout: tuple[float, float] = DEFAULT_TIMEOUT,
) -> Path:
    """Download ``url`` to ``output_path`` atomically, verifying the body.

    The body is streamed to ``<output_path>.part`` and renamed into place
    only after it passed the checks, so ``output_path`` never holds a
    partial file: an existing ``output_path`` can be trusted and skipped by
    the callers, while a leftover ``.part`` is always re-downloaded.

    Parameters
    ----------
    url : str
        URL to fetch.
    output_path : Path
        Final location of the file.
    session : requests.Session or None, optional
        Session with retry logic. If None, a new session is created.
    timeout : tuple[float, float], optional
        ``(connect, read)`` timeout in seconds, by default `DEFAULT_TIMEOUT`.

    Returns
    -------
    Path
        ``output_path``.

    Raises
    ------
    requests.HTTPError
        On a non-2xx status (after the session's retries).
    requests.Timeout or requests.ConnectionError
        When the server does not answer within ``timeout`` (the retry
        adapter reports an exhausted read timeout as a ``ConnectionError``
        whose message says "Read timed out").
    DownloadError
        When the body is an HTML page instead of the file, or shorter than
        the advertised ``Content-Length``.

    """
    output_path = Path(output_path)
    part_path = output_path.with_name(output_path.name + PART_SUFFIX)
    if session is None:
        session = create_session()

    try:
        with session.get(url, stream=True, timeout=timeout) as response:
            response.raise_for_status()
            content_type = response.headers.get("Content-Type", "")
            if "text/html" in content_type.lower():
                msg = f"{url} returned an HTML page ({content_type}), not a data file"
                raise DownloadError(msg)

            expected = _expected_length(response)
            received = 0
            with open(part_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                    if received == 0 and _looks_like_html(chunk):
                        msg = f"{url} returned an HTML page, not a data file"
                        raise DownloadError(msg)
                    f.write(chunk)
                    received += len(chunk)

        if expected is not None and received != expected:
            msg = (
                f"{url}: received {received} bytes but Content-Length is"
                f" {expected}; discarding the partial file"
            )
            raise DownloadError(msg)
    except BaseException:
        part_path.unlink(missing_ok=True)
        raise

    part_path.replace(output_path)
    return output_path


def _expected_length(response: requests.Response) -> int | None:
    """``Content-Length`` of a response whose body is not transfer-encoded."""
    if response.headers.get("Content-Encoding"):
        # iter_content yields decoded bytes; the header counts encoded ones
        return None
    value = response.headers.get("Content-Length")
    try:
        return int(value) if value is not None else None
    except ValueError:
        return None


def _is_staged(output_path: Path) -> bool:
    """Whether a previous run left a complete file at ``output_path``.

    Files are renamed into place only once complete (`download_file`), so an
    existing, non-empty file is complete. A leftover ``.part`` is not.
    """
    part_path = output_path.with_name(output_path.name + PART_SUFFIX)
    if part_path.exists():
        logger.debug(f"Removing incomplete download {part_path}")
        part_path.unlink()
    return output_path.exists() and output_path.stat().st_size > 0


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
    version : {"0.1", "0.2", "0.3"}, optional
        UNR data version. Default is "0.3".
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
    DownloadError
        If the server answered with an HTML page or a truncated body.

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

    # Skip if already downloaded (a leftover .part is not a download)
    if _is_staged(output_path):
        logger.debug(f"Lookup table already exists: {output_path}")
        return output_path

    url = LOOKUP_URL.format(version=version)
    logger.info(f"Downloading lookup table from {url}")

    download_file(url, output_path, session=session)
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
    grid_type: GridType = DEFAULT_GRID_TYPE,
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
    version : {"0.1", "0.2", "0.3"}, optional
        UNR data version. Default is "0.3".
    grid_type : {"constant", "variable"}, optional
        ``"constant"``: precomputed linear rates (Version 0.3, IGS20/NA only).
        ``"variable"``: time-variable positions. Default is "constant".
    session : requests.Session or None, optional
        Session with retry logic. If None, a new session is created.

    Returns
    -------
    Path
        Path to the downloaded file, named by `grid_file_name`.

    Raises
    ------
    ValueError
        If version, or the grid type for this version/plate, is not supported.
    requests.HTTPError
        If download fails.
    DownloadError
        If the server answered with an HTML page or a truncated body.

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

    # Handle IGS14/IGS20 plate compatibility (not available in v0.2+)
    if plate == "IGS14" and version in ("0.2", "0.3"):
        plate = "IGS20"

    if grid_type == "constant" and (
        version not in CONSTANT_VERSIONS or plate not in CONSTANT_PLATES
    ):
        msg = (
            "UNR publishes the constant grid only for versions"
            f" {sorted(CONSTANT_VERSIONS)} and plates {sorted(CONSTANT_PLATES)},"
            f" got version='{version}', plate='{plate}'"
        )
        raise ValueError(msg)

    # Build URL and output path
    filename = f"{plate}/{grid_id:06d}_{plate}.tenv8"
    url = f"{GRID_BASE_URLS[grid_type].format(version=version)}/{filename}"

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / grid_file_name(grid_id, plate, grid_type)

    # Skip if already downloaded (a leftover .part is not a download)
    if _is_staged(output_path):
        return output_path

    return download_file(url, output_path, session=session)


def download_grid_files(
    grid_ids: Sequence[int],
    output_dir: Path,
    plate: PlateType = "IGS20",
    version: VersionType = DEFAULT_VERSION,
    grid_type: GridType = DEFAULT_GRID_TYPE,
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
    version : {"0.1", "0.2", "0.3"}, optional
        UNR data version. Default is "0.3".
    grid_type : {"constant", "variable"}, optional
        UNR grid product. Default is "constant".
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
        grid_type=grid_type,
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


def get_frame_grid_points(frame_gdf, grid_gdf, margin_deg=0.5):
    """Get grid points within a buffered frame boundary using proper projections.

    Parameters
    ----------
    frame_gdf : GeoDataFrame
        GeoDataFrame containing OPERA DISP frame geometry
    grid_gdf : GeoDataFrame
        GeoDataFrame containing grid points
    margin_deg : float, default=0.5
        Buffer margin in degrees (converted to meters based on latitude)

    Returns
    -------
    tuple
        (grid_ids, grid_gdf_filtered) - list of grid IDs and filtered GeoDataFrame

    """
    # Get rough latitude from bounds to determine projection
    bounds = frame_gdf.to_crs(epsg=4326).total_bounds  # [minx, miny, maxx, maxy]
    center_lat = (bounds[1] + bounds[3]) / 2  # Average of min and max latitude

    # Choose appropriate projection based on latitude
    if center_lat > 50:  # Alaska/high latitudes
        target_crs = "EPSG:3338"  # Alaska Albers Equal Area
        margin_m = margin_deg * 111000  # Convert degrees to meters (~111km per degree)
    else:
        # For lower latitudes, use appropriate UTM zone
        target_crs = frame_gdf.estimate_utm_crs()
        margin_m = margin_deg * 111000

    # Reproject to equal-area CRS (avoids distortion issues)
    frame_projected = frame_gdf.to_crs(target_crs)
    grid_projected = grid_gdf.to_crs(target_crs)

    # Buffer the frame polygon in projected coordinates
    buffered = frame_projected.buffer(margin_m)

    # Filter grid points to buffered frame area
    mask = grid_projected.intersects(buffered.iloc[0])
    grid_gdf_filtered = grid_gdf[mask]  # Return in original CRS
    grid_ids = grid_gdf_filtered.index.tolist()

    return grid_ids, grid_gdf_filtered


def download_unr_grid(
    frame_id: int,
    output_dir: Path,
    margin_deg: float = 0.5,
    plate: PlateType = "IGS20",
    version: VersionType = DEFAULT_VERSION,
    grid_type: GridType = DEFAULT_GRID_TYPE,
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
    version : {"0.1", "0.2", "0.3"}, optional
        UNR grid version. Default is "0.3".
    grid_type : {"constant", "variable"}, optional
        UNR grid product. Default is "constant".
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
    ...     lookup_table=data_dir / "grid_latlon_lookup_v0.3.txt",
    ...     data_dir=data_dir,
    ...     frame_id=8882
    ... )
    >>> df = grid.to_dataframe()

    """
    logger.info(
        f"Downloading UNR grid for frame {frame_id} "
        f"(plate={plate}, version={version}, grid_type={grid_type},"
        f" margin={margin_deg}°)"
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
    normalize_longitude = version in ("0.2", "0.3")

    lookup = load_lookup_table(lookup_file, normalize_longitude=normalize_longitude)
    grid_gdf = gpd.GeoDataFrame(
        lookup,
        geometry=gpd.points_from_xy(x=lookup.lon, y=lookup.lat),
        crs="EPSG:4326",
    )

    # Get frame geometry
    frame_gdf = get_frame_geojson([frame_id], as_geodataframe=True)

    # Filter grid points to expanded frame bounds
    grid_ids, _ = get_frame_grid_points(frame_gdf, grid_gdf, margin_deg)

    logger.info(f"Found {len(grid_ids)} grid points within frame bounds")

    # Download .tenv8 files for all grid points in parallel
    download_grid_files(
        grid_ids=grid_ids,
        output_dir=output_dir,
        plate=plate,
        version=version,
        grid_type=grid_type,
        max_workers=max_workers,
    )

    logger.info(f"Download complete. Data saved to {output_dir}")

    return output_dir
