import logging
from datetime import datetime
from pathlib import Path
from typing import Literal

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from geepers.gps_sources import UnrGridSource
from opera_utils import get_frame_geojson

logger = logging.getLogger(__name__)


def download_unr_grid(
    frame_id: int,
    output_dir: Path,
    start: datetime | None = None,
    end: datetime | None = None,
    margin_deg: float = 0.5,
    plate: Literal["NA", "PA", "IGS14", "IGS20"] = "IGS20",
    version: Literal["0.1", "0.2"] = "0.2",
) -> None:
    """Download UNR gridded GNSS timeseries for a given frame.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    output_dir : Path
        Output directory for downloaded data.
    start : datetime or None, optional
        Start date for timeseries. If None, downloads from beginning.
    end : datetime or None, optional
        End date for timeseries. If None, downloads until present.
    margin_deg : float, default=0.5
        Margin in degrees to expand frame bounding box.
    plate : {"NA", "PA", "IGS14", "IGS20"}, default="IGS20"
        Reference plate for velocity computation.
    version : {"0.1", "0.2"}, default="0.2"
        UNR grid version to download.

    """
    logger.info(
        f"Downloading UNR gridded GNSS timeseries for frame {frame_id} "
        f"from {start or 'beginning'} to {end or 'present'}"
    )

    selected_frame = get_frame_geojson([frame_id], as_geodataframe=True)
    frame_bounds = tuple(selected_frame.bounds.values[0])

    extended_bbox = (
        frame_bounds[0] - margin_deg,
        frame_bounds[1] - margin_deg,
        frame_bounds[2] + margin_deg,
        frame_bounds[3] + margin_deg,
    )

    grid = UnrGridSource(version=version)
    ts_grid_df = grid.timeseries_many(bbox=extended_bbox)
    ts_grid_df["date"] = pd.to_datetime(ts_grid_df["date"])

    output_dir.mkdir(parents=True, exist_ok=True)

    _save_parquet(ts_grid_df, frame_id, plate, output_dir)


def _save_parquet(
    df: pd.DataFrame,
    frame_id: int,
    plate: str,
    output_dir: Path,
) -> None:
    """Save timeseries DataFrame to parquet with metadata.

    Parameters
    ----------
    df : pd.DataFrame
        Timeseries data to save.
    frame_id : int
        Frame identifier for metadata.
    plate : str
        Reference plate for metadata.
    output_dir : Path
        Output directory.

    """
    df_no_geom = df.drop(columns=["geometry"], errors="ignore")
    table = pa.Table.from_pandas(df_no_geom)

    metadata = {
        "frame_id": str(frame_id),
        "source": "UNR grid",
        "description": "Time series of east/north/up displacements and uncertainties",
        "reference_frame": plate,
    }

    existing_metadata = table.schema.metadata or {}
    encoded_metadata = {k: v.encode() for k, v in metadata.items()}
    table = table.replace_schema_metadata({**existing_metadata, **encoded_metadata})

    output_file = output_dir / f"unr_grid_frame{frame_id}.parquet"
    pq.write_table(table, output_file)
    logger.info(f"Saved parquet to {output_file}")
