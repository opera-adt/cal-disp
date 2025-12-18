from __future__ import annotations

import logging
from datetime import datetime, timedelta
from pathlib import Path

from opera_utils.disp._download import run_download

# Constants
DEFAULT_NUM_WORKERS = 2
DATE_FORMATS = ("%Y-%m-%d", "%Y%m%d")
SINGLE_DATE_BUFFER_DAYS = 1

logger = logging.getLogger(__name__)


def _adjust_single_date_range(
    start: datetime | None,
    end: datetime | None,
) -> tuple[datetime | None, datetime | None]:
    """Expand date range if start and end are identical.

    When querying for a single date, expand the range by adding buffer days
    to ensure products are captured.

    Parameters
    ----------
    start : datetime or None
        Start datetime.
    end : datetime or None
        End datetime.

    Returns
    -------
    tuple[datetime or None, datetime or None]
        Adjusted start and end datetimes.

    """
    if start and end and start == end:
        logger.info(
            "Single date query detected. Expanding range by "
            f"±{SINGLE_DATE_BUFFER_DAYS} day(s)"
        )
        buffer = timedelta(days=SINGLE_DATE_BUFFER_DAYS)
        return start - buffer, end + buffer
    return start, end


# NOTE expose url type to enable using s3 paths
def download_disp(
    frame_id: int,
    output_dir: Path,
    start: datetime | None = None,
    end: datetime | None = None,
    num_workers: int = DEFAULT_NUM_WORKERS,
) -> None:
    """Download DISP-S1 products for a frame.

    Downloads displacement products from the OPERA DISP-S1 archive for
    the specified frame and date range. Products are filtered based on
    the secondary date of each interferogram.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    output_dir : Path
        Directory where products will be saved.
    start : datetime or None, optional
        Start date for query (based on secondary date).
        Default is None (no start limit).
    end : datetime or None, optional
        End date for query (based on secondary date).
        Default is None (no end limit).
    num_workers : int, optional
        Number of parallel download workers. Default is 2.

    Raises
    ------
    ValueError
        If frame ID is not in the database or no products are found.

    Notes
    -----
    Date queries are based on the secondary (later) date of each
    interferometric pair. If start and end dates are identical,
    the range is automatically expanded by ±1 day and num_workers
    is set to 1 to ensure the specific product is captured.

    Examples
    --------
    >>> download_frame_products(
    ...     frame_id=8887,
    ...     output_dir=Path("./data"),
    ...     start=datetime(2024, 1, 1),
    ...     end=datetime(2024, 12, 31)
    ... )

    """
    logger.info(f"Downloading DISP-S1 products for frame {frame_id}")

    start_adjusted, end_adjusted = _adjust_single_date_range(start, end)
    logger.info(f"Search window: {start_adjusted} - {end_adjusted}")

    workers = 1 if (start and end and start == end) else num_workers

    run_download(
        frame_id,
        start_datetime=start_adjusted,
        end_datetime=end_adjusted,
        output_dir=Path(output_dir),
        num_workers=workers,
    )

    logger.info(f"Download complete: files saved to {output_dir}")
