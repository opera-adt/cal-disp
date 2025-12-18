from __future__ import annotations

import logging
import re
import shutil
from datetime import date, datetime, timedelta
from pathlib import Path

import asf_search as asf
import geopandas as gpd
import numpy as np
import opera_utils
import pandas as pd
import pyproj
import xarray as xr
from scipy.spatial import ConvexHull
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union
from tqdm import tqdm

logger = logging.getLogger(__name__)

# Constants
CSLC_COLLECTION = "C2777443834-ASF"


def get_cslc_amplitude(ds: xr.Dataset) -> xr.DataArray:
    """Compute amplitude from complex VV polarization data.

    Parameters
    ----------
    ds : xr.Dataset
        CSLC dataset containing VV polarization data.

    Returns
    -------
    xr.DataArray
        Amplitude computed as sqrt(real^2 + imag^2).

    """
    vv = ds["VV"]
    real = vv.data.real
    imag = vv.data.imag
    amp = np.sqrt(real**2 + imag**2)
    return xr.DataArray(amp, coords=vv.coords, dims=vv.dims, name="VV_amplitude")


def get_data_bounds(ds: xr.DataArray) -> Polygon:
    """Extract convex hull of valid (non-NaN) data points.

    Parameters
    ----------
    ds : xr.DataArray
        Data array with x_coordinates and y_coordinates.

    Returns
    -------
    Polygon
        Convex hull polygon around valid data points.

    """
    valid_mask = np.isnan(ds)
    iy, ix = np.where(~valid_mask)
    x_clean = ds.x_coordinates[ix].values
    y_clean = ds.y_coordinates[iy].values
    points = np.column_stack([x_clean, y_clean])
    hull = ConvexHull(points)
    hull_points = points[hull.vertices]
    return Polygon(hull_points)


def get_crs(ds: xr.Dataset) -> pyproj.CRS:
    """Extract CRS from CSLC dataset projection attributes.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with projection variable containing spatial_ref attribute.

    Returns
    -------
    pyproj.CRS
        Coordinate reference system.

    """
    return pyproj.CRS.from_wkt(ds.projection.attrs["spatial_ref"])


def extract_burst_info(
    filename: str | Path,
) -> tuple[str | None, int | None, int | None]:
    """Parse track, burst_id, and swath number from CSLC filename.

    Expects format: T{track}_{burst_id}_{swath} where swath is IW1/IW2/IW3.

    Parameters
    ----------
    filename : str or Path
        CSLC filename to parse.

    Returns
    -------
    tuple[str | None, int | None, int | None]
        Track number, burst ID, and swath number (1-3).
        Returns (None, None, None) if pattern not found.

    """
    pattern = r"T(\d+)[-_](\d+)[-_](IW\d)"
    match = re.search(pattern, str(filename))
    if match:
        track = match.group(1)
        burst_id = int(match.group(2))
        swath = match.group(3)
        swath_num = int(swath[-1])
        return track, burst_id, swath_num
    return None, None, None


def search_cslc_bursts(
    burst_ids: list[str],
    sensing_time: datetime,
    time_window_hours: float = 2.0,
) -> asf.ASFSearchResults:
    """Search for CSLC products matching burst IDs and sensing time.

    Parameters
    ----------
    burst_ids : list[str]
        List of OPERA burst IDs (e.g., ['T087_185678_IW1']).
    sensing_time : datetime
        Target sensing time for the acquisition.
    time_window_hours : float, optional
        Time window in hours around sensing_time. Default is 2.0.

    Returns
    -------
    asf.ASFSearchResults
        Search results containing matching CSLC products.

    """
    dt = timedelta(hours=time_window_hours)
    logger.info(
        f"Searching CSLC: time={sensing_time.isoformat()}, "
        f"window=±{time_window_hours}h, bursts={len(burst_ids)}"
    )
    return asf.search(
        operaBurstID=burst_ids,
        start=sensing_time - dt,
        end=sensing_time + dt,
        collections=[CSLC_COLLECTION],
    )


def download_cslc_files(
    results: asf.ASFSearchResults,
    output_dir: Path,
    target_date: date,
    n_processes: int = 5,
) -> list[Path]:
    """Download CSLC files for a specific date.

    Parameters
    ----------
    results : asf.ASFSearchResults
        ASF search results to filter and download.
    output_dir : Path
        Directory to save downloaded files.
    target_date : date
        Date to filter results by startTime.
    n_processes : int, optional
        Number of parallel download processes. Default is 5.

    Returns
    -------
    list[Path]
        Paths to downloaded CSLC files.

    Raises
    ------
    ValueError
        If no results are found for the target date.

    """
    # Convert results to GeoDataFrame and filter by date
    cslc_df = gpd.GeoDataFrame.from_features(results.geojson())
    logger.debug(f"Initial search returned {len(cslc_df)} products")

    if len(cslc_df) == 0:
        raise ValueError("No CSLC products found in search results")

    cslc_df["startTime"] = pd.to_datetime(cslc_df["startTime"])
    cslc_df["date"] = cslc_df["startTime"].dt.date

    cslc_df_date = cslc_df.groupby("date").agg(
        frame_urls=("url", lambda x: list(x)),
        cslc_list=("fileName", lambda x: list(x)),
        n_urls=("url", "count"),
    )

    if target_date not in cslc_df_date.index:
        available_dates = sorted(cslc_df_date.index)
        raise ValueError(
            f"No products found for {target_date}. "
            f"Available dates: {', '.join(str(d) for d in available_dates)}"
        )

    specific_date = cslc_df_date.loc[target_date]
    cslc_df = cslc_df[cslc_df.fileName.isin(specific_date.cslc_list)]
    subset_results = asf.ASFSearchResults([results.data[i] for i in cslc_df.index])

    logger.info(f"Downloading {len(subset_results)} CSLC files for {target_date}")

    # Download to temp directory
    temp_dir = output_dir / "tmp"
    temp_dir.mkdir(exist_ok=True)
    subset_results.download(path=temp_dir, processes=n_processes)

    return [
        temp_dir / subset_results[i].properties["fileName"]
        for i in range(len(subset_results))
    ]


def process_cslc_bounds(cslc_files: list[Path], epsg: int) -> gpd.GeoDataFrame:
    """Extract valid data bounds from CSLC files and reproject to target CRS.

    Parameters
    ----------
    cslc_files : list[Path]
        Paths to CSLC NetCDF files.
    epsg : int
        Target EPSG code for output CRS (e.g., 4326 for WGS84).

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with burst info and geometry bounds in target CRS.

    Notes
    -----
    Each CSLC file may have different native CRS (e.g., different UTM zones).
    All geometries are reprojected to the target EPSG before combining.

    """
    burst_gdfs: list[gpd.GeoDataFrame] = []
    crs_counts: dict[str, int] = {}

    logger.info(f"Processing bounds for {len(cslc_files)} CSLC files")
    for cslc in tqdm(
        cslc_files,
        desc="Processing CSLC bounds",
        disable=logger.level > logging.INFO,
    ):
        # Get amplitude and bounds
        ds = xr.open_dataset(cslc, group="data")
        amp_ds = get_cslc_amplitude(ds)
        bounds = get_data_bounds(amp_ds)

        # Extract burst metadata
        track, burst_id, swath_num = extract_burst_info(cslc.name)

        # Get native CRS for this burst
        native_crs = get_crs(ds)
        crs_str = str(native_crs)
        crs_counts[crs_str] = crs_counts.get(crs_str, 0) + 1

        # Create GeoDataFrame in native CRS
        gdf = gpd.GeoDataFrame(
            [
                {
                    "filename": cslc.name,
                    "track": track,
                    "burst_id": burst_id,
                    "swath": swath_num,
                }
            ],
            geometry=[bounds],
            crs=native_crs,
        )

        # Reproject to target CRS
        gdf = gdf.to_crs(epsg=epsg)
        burst_gdfs.append(gdf)

    # Log CRS diversity
    if len(crs_counts) > 1:
        logger.info(f"Found {len(crs_counts)} different CRS across bursts")
        for crs, count in crs_counts.items():
            logger.debug(f"  {crs}: {count} bursts")
    else:
        logger.debug(f"All bursts share CRS: {list(crs_counts.keys())[0]}")

    logger.debug(f"Reprojected all bursts to EPSG:{epsg}")
    return pd.concat(burst_gdfs, ignore_index=True)


def create_nonoverlapping_tiles(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Create non-overlapping burst tiles with swath/burst priority.

    Removes overlaps with priority: IW1 > IW2 > IW3, and within each swath,
    lower burst_id takes priority over higher burst_id.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame with burst_id, swath, and geometry columns.

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with non-overlapping burst geometries.

    """
    gdf_sorted = gdf.sort_values(["burst_id", "swath"])
    logger.info(f"Creating non-overlapping tiles for {len(gdf_sorted)} bursts")

    burst_polygons: list[Polygon | MultiPolygon] = []
    # Track cumulative coverage per swath
    swath_coverage: dict[int, Polygon | MultiPolygon] = {}

    for swath_num in [1, 2, 3]:
        swath_bursts = gdf_sorted[gdf_sorted["swath"] == swath_num]

        if len(swath_bursts) == 0:
            logger.debug(f"No bursts found for swath IW{swath_num}")
            continue

        swath_bursts = swath_bursts.sort_values("burst_id")
        cumulative_within_swath = None
        logger.debug(f"Processing {len(swath_bursts)} bursts for swath IW{swath_num}")

        for idx, row in swath_bursts.iterrows():
            burst_geom = row.geometry

            # Remove overlaps with all previous swaths
            for prev_swath_num in range(1, swath_num):
                if prev_swath_num in swath_coverage:
                    burst_geom = burst_geom.difference(swath_coverage[prev_swath_num])

            # Remove overlaps with previous bursts in same swath
            if cumulative_within_swath is not None:
                burst_geom = burst_geom.difference(cumulative_within_swath)

            burst_polygons.append(
                {
                    "swath": row.swath,
                    "burst_id": row.burst_id,
                    "filename": row.filename,
                    "geometry": burst_geom,
                }
            )

            # Update cumulative coverage
            if cumulative_within_swath is None:
                cumulative_within_swath = row.geometry
            else:
                cumulative_within_swath = unary_union(
                    [cumulative_within_swath, row.geometry]
                )

        swath_coverage[swath_num] = cumulative_within_swath

    return gpd.GeoDataFrame(burst_polygons, crs=gdf.crs)


def generate_s1_burst_tiles(
    frame_id: int,
    sensing_time: datetime,
    output_dir: Path,
    time_window_hours: float = 2.0,
    n_download_processes: int = 5,
) -> Path:
    """Generate non-overlapping burst tiles for a frame and sensing time.

    Downloads CSLC data, processes bursts to create non-overlapping polygons,
    and saves to GeoJSON. Priority: IW1 > IW2 > IW3, lower burst_id first.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    sensing_time : datetime
        Sensing time to search for CSLC products.
    output_dir : Path
        Directory to save output GeoJSON and temporary files.
    time_window_hours : float, optional
        Time window in hours for searching CSLC products. Default is 2.0.
    n_download_processes : int, optional
        Number of parallel download processes. Default is 5.

    Returns
    -------
    Path
        Path to the generated GeoJSON file.

    Raises
    ------
    ValueError
        If no bursts are found or download fails.

    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    logger.info(f"Starting burst tile generation for frame {frame_id}")

    # Get burst IDs and EPSG for frame
    burst_ids = opera_utils.get_burst_ids_for_frame(frame_id)
    burst_ids = [b.upper() for b in burst_ids]
    logger.info(f"Found {len(burst_ids)} burst IDs for frame {frame_id}")

    epsg = opera_utils.get_frame_bbox(frame_id)[0]
    logger.debug(f"Target CRS for frame: EPSG:{epsg}")

    # Search for CSLC products
    results = search_cslc_bursts(burst_ids, sensing_time, time_window_hours)
    logger.info(f"Found {len(results)} CSLC products")

    # Download files
    target_date = sensing_time.date()
    cslc_files = download_cslc_files(
        results, output_dir, target_date, n_download_processes
    )
    logger.info(f"Downloaded {len(cslc_files)} CSLC files")

    # Process bounds
    cslc_gdf = process_cslc_bounds(cslc_files, epsg=epsg)

    # Create non-overlapping tiles
    burst_gdf = create_nonoverlapping_tiles(cslc_gdf)

    # Log summary statistics
    swath_counts = burst_gdf.groupby("swath").size()
    logger.info(f"Generated {len(burst_gdf)} non-overlapping tiles")
    for swath_num, count in swath_counts.items():
        logger.info(f"  IW{swath_num}: {count} tiles")

    # Save to GeoJSON
    output_file = output_dir / f"{target_date}_tiles.geojson"
    burst_gdf.to_file(output_file)
    logger.info(f"Saved tiles to {output_file}")

    # Clean up temp directory
    temp_dir = output_dir / "tmp"
    if temp_dir.exists():
        shutil.rmtree(temp_dir)
        logger.debug(f"Cleaned up temporary directory: {temp_dir}")

    return output_file
