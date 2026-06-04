"""Generate an event mask GeoTIFF from a GeoJSON and an OPERA DISP-S1 product.

Input GeoJSON
-------------
Each feature must have the following fields:

    - ``id``         : unique feature identifier
    - ``frame_id``   : OPERA frame identifier (integer)
    - ``event_date`` : ISO-8601 date string (YYYY-MM-DD or YYYY-MM-DDTHH:MM:SS)
    - ``geometry``   : polygon geometry in the file's declared CRS

Feature selection
-----------------
The OPERA DISP-S1 product filename encodes the frame ID and acquisition date
range (date1, date2).  A GeoJSON feature is selected when:

    - its ``frame_id`` matches the product frame ID, AND
    - its ``event_date`` falls within the product epoch [date1, date2]

Coordinate reference system
----------------------------
- Input GeoJSON    : any CRS (auto-detected; EPSG:4326 assumed if undeclared)
- Output GeoTIFF   : same CRS as the OPERA DISP-S1 product (variable UTM zone)

Geometries are reprojected from the GeoJSON CRS to the product CRS before
rasterization, so the mask is pixel-for-pixel aligned with the input product.

Output convention
-----------------
- ``1`` : valid pixel (no event)
- ``0`` : masked pixel (event region)

Output naming
-------------
The output filename is derived automatically from both inputs:

    <product_stem>_<geojson_stem>_mask.tif

Examples:
    events.geojson          -> ..._events_mask.tif
    continuous_defo.geojson -> ..._continuous_defo_mask.tif

Use ``--output`` to override this path.

Usage
-----
::

    # events mask — output auto-named from GeoJSON stem
    python generate_event_mask.py events.geojson \
        OPERA_L3_DISP-S1_IW_F36540_VV_20160724T015809Z_20160805T015809Z_v1.0.nc

    # continuous deformation mask
    python generate_event_mask.py continuous_defo.geojson \
        OPERA_L3_DISP-S1_IW_F36540_VV_20160724T015809Z_20160805T015809Z_v1.0.nc

    # explicit output path
    python generate_event_mask.py events.geojson product.nc --output masks/my_mask.tif

"""

from __future__ import annotations

import argparse
import logging
import re
from datetime import date, datetime
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger(__name__)

# Matches _F<digits>_ anywhere in the filename
_FRAME_RE = re.compile(r"_F(\d+)_")

# Matches the first two ISO date-time stamps:  _<YYYYMMDDTHHMMSSZ>_<YYYYMMDDTHHMMSSZ>_
_DATES_RE = re.compile(r"_(\d{8}T\d{6}Z)_(\d{8}T\d{6}Z)_")

_DATE_FMT = "%Y%m%dT%H%M%SZ"


def parse_product_metadata(product_path: Path) -> tuple[int, date, date]:
    """Parse frame ID and epoch date range from an OPERA DISP-S1 filename.

    Parameters
    ----------
    product_path : Path
        Path to the OPERA DISP-S1 NetCDF file.

    Returns
    -------
    tuple
        ``(frame_id, date1, date2)``

    Raises
    ------
    ValueError
        If the filename does not contain a recognisable frame ID or date pair.

    Examples
    --------
    ::

        frame_id, d1, d2 = parse_product_metadata(
            Path("OPERA_L3_DISP-S1_IW_F36540_VV_20160724T015809Z_20160805T015809Z_v1.0.nc")
        )
        # frame_id=36540, d1=date(2016,7,24), d2=date(2016,8,5)

    """
    name = product_path.name

    frame_match = _FRAME_RE.search(name)
    if not frame_match:
        msg = f"Cannot parse frame ID from filename: {name!r}"
        raise ValueError(msg)
    frame_id = int(frame_match.group(1))

    date_match = _DATES_RE.search(name)
    if not date_match:
        msg = f"Cannot parse date pair from filename: {name!r}"
        raise ValueError(msg)
    date1 = datetime.strptime(date_match.group(1), _DATE_FMT).date()
    date2 = datetime.strptime(date_match.group(2), _DATE_FMT).date()

    return frame_id, date1, date2


def _read_product_georef(product_path: Path) -> tuple[Any, Any, int, int]:
    """Read CRS, affine transform, width, and height from a DISP-S1 NetCDF.

    Parameters
    ----------
    product_path : Path
        Path to the OPERA DISP-S1 NetCDF file.

    Returns
    -------
    tuple
        ``(crs_wkt, transform, width, height)``

    """
    import xarray as xr
    from rasterio.transform import from_bounds

    with xr.open_dataset(product_path) as ds:
        # CRS — prefer crs_wkt from spatial_ref variable
        crs = "EPSG:4326"
        if "spatial_ref" in ds:
            sr_attrs = ds["spatial_ref"].attrs
            for attr in ("crs_wkt", "spatial_ref", "wkt"):
                if attr in sr_attrs:
                    crs = sr_attrs[attr]
                    break
        elif "crs" in ds.attrs:
            crs = ds.attrs["crs"]
        else:
            logger.warning("No CRS found in product, assuming EPSG:4326")

        x = ds.coords["x"].values
        y = ds.coords["y"].values

    width = len(x)
    height = len(y)
    x_res = abs(x[1] - x[0])
    y_res = abs(y[1] - y[0])

    transform = from_bounds(
        x.min() - x_res / 2,
        y.min() - y_res / 2,
        x.max() + x_res / 2,
        y.max() + y_res / 2,
        width,
        height,
    )

    return crs, transform, width, height


def default_output_path(product_path: Path, geojson_path: Path) -> Path:
    """Return the default output mask path derived from the product and GeoJSON paths.

    The output is placed alongside the product and named
    ``<product_stem>_<geojson_stem>_mask.tif``, so the GeoJSON filename
    determines the mask label automatically.

    Parameters
    ----------
    product_path : Path
        Path to the OPERA DISP-S1 NetCDF product.
    geojson_path : Path
        Path to the input GeoJSON file.

    Returns
    -------
    Path
        Derived output path, e.g. for ``events.geojson``:
        ``OPERA_L3_DISP-S1_IW_F36540_..._events_mask.tif``

        For ``continuous_defo.geojson``:
        ``OPERA_L3_DISP-S1_IW_F36540_..._continuous_defo_mask.tif``

    """
    return product_path.with_name(f"{product_path.stem}_{geojson_path.stem}_mask.tif")


def generate_event_mask(
    product_path: Path,
    geojson_path: Path,
    output_path: Path | None = None,
    filter_by_date: bool = True,
) -> Path:
    """Generate an event mask GeoTIFF aligned to an OPERA DISP-S1 product.

    Reads event polygons from ``geojson_path``, selects those whose ``frame_id``
    and ``event_date`` match the product epoch, reprojects them to the product
    CRS, and rasterizes them onto a grid that is pixel-for-pixel identical to
    the input product.

    The output filename is derived automatically from both input filenames as
    ``<product_stem>_<geojson_stem>_mask.tif`` unless ``output_path`` is given
    explicitly.

    Parameters
    ----------
    product_path : Path
        Path to the OPERA DISP-S1 NetCDF product.
    geojson_path : Path
        Path to the GeoJSON file containing event or deformation polygons.
        Required feature properties: ``frame_id`` (int), ``geometry``.
        ``event_date`` is additionally required when ``filter_by_date=True``.
    output_path : Path, optional
        Destination path for the output GeoTIFF mask.  When omitted the file is
        written next to the product as ``<product_stem>_<geojson_stem>_mask.tif``.
    filter_by_date : bool, optional
        When ``True`` (default) only features whose ``event_date`` falls within
        the product epoch ``[date1, date2]`` are masked.  Set to ``False`` for
        continuous deformation databases where all features matching the
        ``frame_id`` should always be masked regardless of date.

    Returns
    -------
    Path
        Path to the written GeoTIFF mask.

    Examples
    --------
    ::

        out = generate_event_mask(
            product_path=Path("OPERA_L3_DISP-S1_IW_F36540_VV_20160724T015809Z_20160805T015809Z_v1.0.nc"),
            geojson_path=Path("events.geojson"),
        )
        print(out)
        # OPERA_L3_DISP-S1_IW_F36540_VV_20160724T015809Z_20160805T015809Z
        # _v1.0_events_mask.tif

        out = generate_event_mask(
            product_path=Path("OPERA_L3_DISP-S1_IW_F36540_VV_20160724T015809Z_20160805T015809Z_v1.0.nc"),
            geojson_path=Path("continuous_defo.geojson"),
        )
        print(out)
        # OPERA_L3_DISP-S1_IW_F36540_VV_20160724T015809Z_20160805T015809Z
        # _v1.0_continuous_defo_mask.tif

    """
    import geopandas as gpd
    import pandas as pd
    import rasterio as rio
    from rasterio.features import rasterize

    if output_path is None:
        output_path = default_output_path(product_path, geojson_path)

    frame_id, date1, date2 = parse_product_metadata(product_path)
    logger.info(f"Product  : frame={frame_id}, epoch={date1} -> {date2}")

    crs, transform, width, height = _read_product_georef(product_path)

    # Load GeoJSON and normalise types.
    # Respect whatever CRS the file declares.  Fall back to EPSG:4326 only when
    # the file omits the CRS member, which geopandas would otherwise leave as None
    # and which RFC 7946 defines as the GeoJSON default.
    gdf = gpd.read_file(geojson_path)
    if gdf.crs is None:
        logger.warning(
            "GeoJSON has no CRS declaration; assuming EPSG:4326 (RFC 7946 default)"
        )
        gdf = gdf.set_crs("EPSG:4326")
    geojson_crs = gdf.crs.to_string()
    logger.info(f"GeoJSON CRS: {geojson_crs}")

    # Validate required columns before touching them
    required = {"frame_id", "event_date"} if filter_by_date else {"frame_id"}
    missing_cols = required - set(gdf.columns)
    if missing_cols:
        raise ValueError(
            f"GeoJSON {geojson_path.name!r} is missing required column(s): "
            f"{sorted(missing_cols)}.  Each feature must have: "
            f"{sorted(required)}."
        )

    gdf["frame_id"] = gdf["frame_id"].astype(int)

    if filter_by_date:
        gdf["event_date"] = pd.to_datetime(gdf["event_date"]).dt.date
        matched = gdf[
            (gdf["frame_id"] == frame_id)
            & (gdf["event_date"] >= date1)
            & (gdf["event_date"] <= date2)
        ]
    else:
        # Continuous deformation: mask for any epoch of this frame
        matched = gdf[gdf["frame_id"] == frame_id]
    logger.info(f"Matched  : {len(matched)} event feature(s)")

    # Start with all-valid mask
    mask = np.ones((height, width), dtype=np.uint8)

    if len(matched) > 0:
        # Reproject from the GeoJSON CRS to the product CRS before rasterizing
        # so the burned pixels align with the product grid exactly.
        matched_proj = matched.to_crs(crs)
        logger.info(f"Reprojected geometries from {geojson_crs} to product CRS")
        shapes = [
            (geom, 0)
            for geom in matched_proj.geometry
            if geom is not None and not geom.is_empty
        ]

        if shapes:
            mask = rasterize(
                shapes,
                out_shape=(height, width),
                transform=transform,
                fill=1,
                dtype=np.uint8,
            )
            logger.info(f"Burned   : {len(shapes)} polygon(s) into mask")
        else:
            logger.warning("All matched geometries were null or empty; mask is all 1")

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with rio.open(
        output_path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype=np.uint8,
        crs=crs,
        transform=transform,
        compress="lzw",
    ) as dst:
        dst.write(mask[np.newaxis, :, :])
        dst.set_band_description(1, "Event mask (1=valid, 0=event region)")

    n_masked = int((mask == 0).sum())
    logger.info(f"Output   : {output_path}  ({n_masked:,} masked pixels)")
    return output_path


def main() -> None:
    """Parse arguments and run the event mask generator."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate an event mask GeoTIFF from a GeoJSON and an OPERA DISP-S1"
            " product.\n\nGeoJSON fields required: id, frame_id, event_date,"
            " geometry.\nAny GeoJSON CRS is accepted; output is always in the product"
            " CRS.\nValues: 1 = valid pixel, 0 = event region."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Output naming (automatic):\n  events.geojson          ->"
            " <product_stem>_events_mask.tif\n  continuous_defo.geojson ->"
            " <product_stem>_continuous_defo_mask.tif\n\nExamples:\n  python"
            " generate_event_mask.py events.geojson product.nc\n  python"
            " generate_event_mask.py continuous_defo.geojson product.nc\n  python"
            " generate_event_mask.py events.geojson product.nc -o masks/custom.tif\n"
        ),
    )
    parser.add_argument(
        "geojson",
        type=Path,
        help=(
            "Path to GeoJSON file with event or deformation polygons (any CRS;"
            " EPSG:4326 assumed if undeclared)."
        ),
    )
    parser.add_argument(
        "product",
        type=Path,
        help="Path to OPERA DISP-S1 NetCDF product.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help=(
            "Override the output GeoTIFF path. By default the mask is written "
            "next to the product as <product_stem>_<geojson_stem>_mask.tif."
        ),
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    generate_event_mask(
        product_path=args.product,
        geojson_path=args.geojson,
        output_path=args.output,  # None → auto-named from product stem
    )


if __name__ == "__main__":
    main()
