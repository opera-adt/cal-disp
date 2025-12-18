from __future__ import annotations

import logging
from datetime import datetime
from pathlib import Path

logger = logging.getLogger(__name__)


def extract_sensing_times_from_file(disp_file: Path) -> list[datetime]:
    """Extract sensing times from DISP-S1 NetCDF filename.

    Parses DISP-S1 filename to extract reference and secondary dates.
    Expected format:
    OPERA_L3_DISP-S1_IW_F{frame}_VV_{ref_date}_{sec_date}_v{version}_{prod_date}.nc

    Parameters
    ----------
    disp_file : Path
        Path to DISP-S1 NetCDF file.

    Returns
    -------
    list[datetime]
        List of unique sensing times from reference and secondary dates.

    Raises
    ------
    FileNotFoundError
        If DISP file does not exist.
    ValueError
        If dates cannot be parsed from filename.

    Examples
    --------
    >>> from pathlib import Path
    >>> filename = Path("OPERA_L3_*.nc")
    >>> times = extract_sensing_times_from_file(filename)  # doctest: +SKIP
    >>> len(times)  # doctest: +SKIP
    2

    """
    if not disp_file.exists():
        msg = f"DISP file not found: {disp_file}"
        raise FileNotFoundError(msg)

    logger.info(f"Extracting sensing times from {disp_file.name}")

    filename = disp_file.name
    parts = filename.split("_")

    # Find parts matching datetime format (YYYYMMDDTHHMMSSZ)
    datetime_parts = [p for p in parts if len(p) == 16 and p.endswith("Z") and "T" in p]

    if len(datetime_parts) < 2:
        msg = (
            f"Cannot parse reference and secondary dates from filename: {filename}. "
            "Expected format: "
            "OPERA_L3_DISP-S1_IW_F{{frame}}_VV_{{ref_date}}_{{sec_date}}_"
            "v{{version}}_{{prod_date}}.nc"
        )
        raise ValueError(msg)

    ref_date_str = datetime_parts[0]
    sec_date_str = datetime_parts[1]

    try:
        ref_date = datetime.strptime(ref_date_str, "%Y%m%dT%H%M%SZ")
        sec_date = datetime.strptime(sec_date_str, "%Y%m%dT%H%M%SZ")
    except ValueError as e:
        msg = f"Failed to parse dates from filename {filename}: {e}"
        raise ValueError(msg) from e

    sensing_times = sorted({ref_date, sec_date})
    logger.info(
        f"Parsed sensing times: {ref_date.isoformat()} (ref), "
        f"{sec_date.isoformat()} (sec)"
    )

    return sensing_times
