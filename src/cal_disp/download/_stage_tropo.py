from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

import asf_search as asf
import geopandas as gpd
import pandas as pd

# Constants
TROPO_COLLECTION = "C3717139408-ASF"

logger = logging.getLogger(__name__)

# Suppress verbose logging from third-party libraries
logging.getLogger("asf_search").setLevel(logging.ERROR)
logging.getLogger("urllib3").setLevel(logging.WARNING)


def find_nearest_scenes(
    target_time: datetime | pd.Timestamp,
    collection: str = TROPO_COLLECTION,
    window_hours: float = 12,
    num_scenes: int = 2,
) -> asf.ASFSearchResults:
    """Find ASF scenes nearest to target time.

    When num_scenes=2, returns one scene before and one after for temporal
    interpolation. When num_scenes=1, returns single nearest scene.

    Parameters
    ----------
    target_time : datetime or pd.Timestamp
        Target sensing time
    collection : str
        ASF collection ID (default: GACOS tropospheric)
    window_hours : float
        Search window ± hours around target
    num_scenes : int
        Number of scenes to return (1 or 2)

    Returns
    -------
    asf_search.ASFSearchResults
        For num_scenes=2: [before, after] bracketing target time.
        For num_scenes=1: single nearest scene.

    """
    t = pd.Timestamp(target_time, tz="UTC")
    dt = pd.Timedelta(hours=window_hours)

    results = asf.search(
        collections=[collection],
        start=t - dt,
        end=t + dt,
    )

    if len(results) == 0:
        return results

    gdf = gpd.GeoDataFrame.from_features(results.geojson())
    gdf["startTime"] = pd.to_datetime(gdf["startTime"]).dt.tz_convert("UTC")
    gdf["result_obj"] = list(results)

    if num_scenes == 1:
        gdf["dt_sec"] = (gdf["startTime"] - t).abs().dt.total_seconds()
        nearest = gdf.loc[gdf["dt_sec"].idxmin()]
        return asf.ASFSearchResults([nearest["result_obj"]])

    # Get one scene before and one after target
    before = gdf[gdf["startTime"] < t]
    after = gdf[gdf["startTime"] >= t]

    selected = []

    if len(before) > 0:
        time_diff = (t - before["startTime"]).dt.total_seconds()
        nearest_before = before.loc[time_diff.idxmin()]
        selected.append(nearest_before["result_obj"])

    if len(after) > 0:
        time_diff = (after["startTime"] - t).dt.total_seconds()
        nearest_after = after.loc[time_diff.idxmin()]
        selected.append(nearest_after["result_obj"])

    # Fall back to N nearest if we can't bracket
    if len(selected) < num_scenes:
        gdf["dt_sec"] = (gdf["startTime"] - t).abs().dt.total_seconds()
        nearest_rows = gdf.nsmallest(num_scenes, "dt_sec")
        selected = nearest_rows["result_obj"].tolist()

    return asf.ASFSearchResults(selected)


def download_tropo(
    disp_times: list[datetime | pd.Timestamp],
    output_dir: Path | str,
    num_workers: int = 4,
    interp: bool = True,
) -> None:
    """Download tropospheric correction data for displacement times.

    Parameters
    ----------
    disp_times : list of datetime or pd.Timestamp
        Displacement measurement times
    output_dir : Path or str
        Output directory for downloads
    num_workers : int
        Parallel download workers
    interp : bool
        If True, get 2 scenes per time (for interpolation).
        If False, get single nearest scene.

    """
    out = Path(output_dir)
    out.mkdir(exist_ok=True, parents=True)

    n_scenes = 2 if interp else 1
    all_scenes = []

    for t in disp_times:
        scenes = find_nearest_scenes(t, num_scenes=n_scenes)
        all_scenes.extend(scenes)

    combined = asf.ASFSearchResults(all_scenes)
    logger.info(f"Downloading {len(combined)} scenes to {out}")
    combined.download(path=out, processes=num_workers)
