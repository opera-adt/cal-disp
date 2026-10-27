"""Build the plate_motion group of calibration products."""

from __future__ import annotations

import xarray as xr

from ._main import set_coord_attrs
from ._utils import make_spatial_ref


def build_plate_motion_dataset(
    plate_motion: xr.DataArray,
    plate: str,
    model: str,
    pair_span_years: float,
    spatial_ref: xr.DataArray | None,
) -> xr.Dataset:
    """Build the ``plate_motion`` dataset: rigid plate motion over the pair in LOS.

    Parameters
    ----------
    plate_motion : xr.DataArray
        LOS displacement (m) of the plate's rigid rotation over the pair, on the
        DISP grid, dimensions (time, y, x).
    plate : str
        Plate of the frame (frame-parameter table), e.g. "NA".
    model : str
        Plate-motion model, e.g. "ITRF2020-PMM".
    pair_span_years : float
        Pair span (secondary - reference, decimal years).
    spatial_ref : xr.DataArray or None
        Spatial reference of the DISP grid.

    Returns
    -------
    xr.Dataset
        Dataset with ``plate_motion`` and ``spatial_ref``.

    """
    da = plate_motion.copy()
    da.attrs.update(
        {
            "units": "meters",
            "long_name": "Plate motion in LOS over the pair",
            "description": (
                f"LOS displacement of rigid {plate} plate rotation ({model}) between"
                " the reference and secondary dates (positive = toward the"
                " satellite, as DISP). DISP - calibration is in IGS20;"
                " DISP - calibration - plate_motion is plate-fixed."
            ),
            "grid_mapping": "spatial_ref",
        }
    )
    data_vars = {"plate_motion": da}
    if spatial_ref is not None:
        data_vars["spatial_ref"] = make_spatial_ref(
            spatial_ref, da.x.values, da.y.values
        )
    ds = set_coord_attrs(xr.Dataset(data_vars))
    ds.attrs.update(
        {
            "group_id": "plate_motion",
            "description": (
                "Rigid plate motion over the pair projected on the line of sight,"
                " to convert the IGS20-calibrated displacement to a plate-fixed"
                " frame: DISP - calibration - plate_motion"
            ),
            "plate": plate,
            "plate_motion_model": model,
            "pair_span_years": float(pair_span_years),
            "units": "meters",
            "Conventions": "CF-1.8",
            "title": "Plate Motion",
        }
    )
    return ds
