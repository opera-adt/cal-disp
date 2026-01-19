"""Utility functions for calibration products."""

from datetime import datetime

import numpy as np
import xarray as xr
from rasterio.transform import Affine

from .._disp import DispProduct


def build_filename(
    disp_product: DispProduct,
    sensor: str,
    version: str,
    production_date: datetime,
) -> str:
    """Build OPERA-compliant filename for calibration product."""
    return (
        f"OPERA_L4_CAL-DISP-{sensor}_"
        f"{disp_product.mode}_"
        f"F{disp_product.frame_id:05d}_"
        f"{disp_product.polarization}_"
        f"{disp_product.reference_date:%Y%m%dT%H%M%S}Z_"
        f"{disp_product.secondary_date:%Y%m%dT%H%M%S}Z_"
        f"v{version}_"
        f"{production_date:%Y%m%dT%H%M%S}Z.nc"
    )


def get_transform(ds: xr.Dataset) -> Affine:
    """Extract affine transform from dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with spatial_ref containing GeoTransform.

    Returns
    -------
    Affine
        Affine transformation.

    Raises
    ------
    ValueError
        If GeoTransform not found.

    """
    gt = ds.spatial_ref.attrs.get("GeoTransform")
    if gt is None:
        raise ValueError("No GeoTransform found in spatial_ref")

    vals = [float(x) for x in gt.split()]
    return Affine(vals[1], vals[2], vals[0], vals[4], vals[5], vals[3])


def get_crs(ds: xr.Dataset) -> str:
    """Extract CRS from dataset.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset with spatial_ref containing crs_wkt.

    Returns
    -------
    str
        CRS as WKT string.

    Raises
    ------
    ValueError
        If crs_wkt not found.

    """
    crs_wkt = ds.spatial_ref.attrs.get("crs_wkt")
    if crs_wkt is None:
        raise ValueError("No crs_wkt found in spatial_ref")
    return crs_wkt


def compute_transform_from_coords(x: np.ndarray, y: np.ndarray) -> Affine:
    """Compute affine transform from coordinate arrays."""
    return Affine.translation(float(x[0]), float(y[0])) * Affine.scale(
        float(x[1] - x[0]),
        float(y[1] - y[0]),
    )


def compute_stats(data: np.ndarray) -> dict[str, float] | None:
    """Compute summary statistics for array, ignoring NaN values."""
    valid_data = data[~np.isnan(data)]

    if len(valid_data) == 0:
        return None

    return {
        "mean": float(np.mean(valid_data)),
        "std": float(np.std(valid_data)),
        "min": float(np.min(valid_data)),
        "max": float(np.max(valid_data)),
    }
