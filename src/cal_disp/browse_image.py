"""Module for creating browse images for the output product."""

from __future__ import annotations

import cmap as cm
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from numpy.typing import ArrayLike
from scipy import ndimage

from cal_disp._types import PathOrStr

DEFAULT_CMAP = cm.Colormap("vik").to_mpl()
DEFAULT_MAX_DIM = 2048
DEFAULT_PERCENTILE = 98.0


def _resize_to_max_pixel_dim(
    arr: ArrayLike, max_dim_allowed: int = DEFAULT_MAX_DIM
) -> np.ndarray:
    """Shrink `arr` so that its *largest* dimension is at most `max_dim_allowed`.

    The aspect ratio is kept; arrays already within the limit are returned
    as a copy (never upsampled). NaNs are masked before the (bilinear) zoom
    and restored afterwards so they do not bleed zeros into the valid data.
    The input array is not modified.
    """
    if max_dim_allowed < 1:
        raise ValueError(f"{max_dim_allowed} is not a valid max image dimension")
    arr = np.asarray(arr, dtype=np.float32)
    scaling_ratio = min(max_dim_allowed / xy for xy in arr.shape)
    if scaling_ratio >= 1:
        return arr.copy()

    nan_mask = np.isnan(arr)
    filled = np.where(nan_mask, 0.0, arr).astype(np.float32)
    valid = (~nan_mask).astype(np.float32)
    # Mask-aware bilinear resampling: mean of the valid input pixels only
    data = ndimage.zoom(filled, scaling_ratio, order=1)
    weight = ndimage.zoom(valid, scaling_ratio, order=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = data / weight
    out[weight < 0.5] = np.nan
    return out.astype(np.float32)


def symmetric_limits(
    arr: ArrayLike, percentile: float = DEFAULT_PERCENTILE
) -> tuple[float, float]:
    """Colour limits ``(-v, v)`` with ``v`` the `percentile` of ``|arr|`` (NaN-free).

    Robust to outliers and centred on zero, so positive and negative
    calibration values get the same visual weight. Falls back to ``(-1, 1)``
    mm if the array has no finite values or is all zero.
    """
    finite = np.abs(np.asarray(arr, dtype=np.float32))
    finite = finite[np.isfinite(finite)]
    v = float(np.percentile(finite, percentile)) if finite.size else 0.0
    if not v > 0:
        v = 1e-3
    return (-v, v)


def _save_to_disk_as_color(
    arr: ArrayLike, fname: PathOrStr, cmap: str, vmin: float, vmax: float
) -> None:
    """Save image array as color to file."""
    plt.imsave(fname, arr, cmap=cmap, vmin=vmin, vmax=vmax)


def make_browse_image_from_arr(
    output_filename: PathOrStr,
    arr: ArrayLike,
    mask: ArrayLike | None = None,
    max_dim_allowed: int = DEFAULT_MAX_DIM,
    cmap: str = DEFAULT_CMAP,
    vmin: float | None = None,
    vmax: float | None = None,
    percentile: float = DEFAULT_PERCENTILE,
) -> tuple[float, float]:
    """Create a PNG browse image for the output product from given array.

    Parameters
    ----------
    output_filename : PathOrStr
        Output PNG.
    arr : ArrayLike
        2-D array (not modified).
    mask : ArrayLike, optional
        Valid-pixel mask (0/False = masked out, shown transparent).
    max_dim_allowed : int
        Largest allowed image dimension (pixels).
    cmap : str
        Matplotlib colormap.
    vmin, vmax : float, optional
        Colour limits. By default symmetric about zero at the `percentile`
        of ``|arr|`` (see ``symmetric_limits``).
    percentile : float
        Percentile used for the default colour limits.

    Returns
    -------
    tuple[float, float]
        The ``(vmin, vmax)`` used.

    """
    arr = np.array(arr, dtype=np.float32, copy=True)
    if mask is not None:
        arr[np.asarray(mask) == 0] = np.nan
    if vmin is None or vmax is None:
        lo, hi = symmetric_limits(arr, percentile)
        vmin = lo if vmin is None else vmin
        vmax = hi if vmax is None else vmax
    arr = _resize_to_max_pixel_dim(arr, max_dim_allowed)
    _save_to_disk_as_color(arr, output_filename, cmap, vmin, vmax)
    return (vmin, vmax)


def make_browse_image_from_nc(
    output_filename: PathOrStr,
    input_filename: PathOrStr,
    max_dim_allowed: int = DEFAULT_MAX_DIM,
    cmap: str = DEFAULT_CMAP,
    vmin: float | None = None,
    vmax: float | None = None,
    percentile: float = DEFAULT_PERCENTILE,
) -> tuple[float, float]:
    """Create a PNG browse image of the ``calibration`` layer of a product."""
    with xr.open_dataset(input_filename, engine="h5netcdf") as ds:
        da = ds["calibration"]
        # Drop time dimension if it exists
        if "time" in da.dims:
            da = da.isel(time=0)
        arr = da.values

    return make_browse_image_from_arr(
        output_filename,
        arr,
        mask=None,
        max_dim_allowed=max_dim_allowed,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        percentile=percentile,
    )
