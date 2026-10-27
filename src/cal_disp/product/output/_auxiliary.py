"""Build auxiliary group for calibration products."""

from __future__ import annotations

import xarray as xr

from ._main import set_coord_attrs
from ._utils import make_spatial_ref


def build_auxiliary_dataset(
    model_3d: dict[str, xr.DataArray] | None,
    model_3d_std: dict[str, xr.DataArray] | None,
    spatial_ref: xr.DataArray | None,
    pair_span_years: float | None = None,
    gnss_sampled: bool = True,
) -> xr.Dataset:
    """Build auxiliary dataset with 3D displacement components.

    Parameters
    ----------
    model_3d : dict[str, xr.DataArray] or None
        3D displacement components (coarse resolution) with keys:
        "north_south", "east_west", "up_down".
    model_3d_std : dict[str, xr.DataArray] or None
        3D displacement uncertainties (coarse resolution).
    spatial_ref : xr.DataArray or None
        Spatial reference data variable (CRS attributes are kept; the
        ``GeoTransform`` is set for the coarse grid).
    pair_span_years : float or None, optional
        Time span of the pair (secondary - reference, decimal years); recorded
        so that velocity = displacement / span.
    gnss_sampled : bool, optional
        False when the components are a zero placeholder (grid types without a
        Venti sampling path); recorded in the group attributes.

    Returns
    -------
    xr.Dataset
        Auxiliary dataset with 3D displacement model.

    """
    data_vars: dict[str, xr.DataArray] = {}

    # Component descriptions
    component_descriptions = {
        "north_south": "GNSS north-south displacement over the pair (positive = north)",
        "east_west": "GNSS east-west displacement over the pair (positive = east)",
        "up_down": "GNSS up-down displacement over the pair (positive = up)",
    }

    std_descriptions = {
        "north_south_std": "1-sigma uncertainty of the north-south displacement",
        "east_west_std": "1-sigma uncertainty of the east-west displacement",
        "up_down_std": "1-sigma uncertainty of the up-down displacement",
    }

    # Add 3D displacement components
    if model_3d:
        for comp in ["north_south", "east_west", "up_down"]:
            if comp in model_3d:
                da = model_3d[comp].copy()
                da.attrs.update(
                    {
                        "description": component_descriptions[comp],
                        "long_name": comp.replace("_", " ").title(),
                        "units": "meters",
                        "grid_mapping": "spatial_ref",
                    }
                )
                data_vars[comp] = da

    # Add 3D displacement uncertainties
    if model_3d_std:
        for comp in ["north_south_std", "east_west_std", "up_down_std"]:
            if comp in model_3d_std:
                da = model_3d_std[comp].copy()
                da.attrs.update(
                    {
                        "description": std_descriptions[comp],
                        "long_name": comp.replace("_", " ").title(),
                        "units": "meters",
                        "grid_mapping": "spatial_ref",
                    }
                )
                data_vars[comp] = da

    # Spatial reference for the coarse grid
    if spatial_ref is not None and data_vars:
        first = next(iter(data_vars.values()))
        data_vars["spatial_ref"] = make_spatial_ref(
            spatial_ref, first.x.values, first.y.values
        )

    ds = set_coord_attrs(xr.Dataset(data_vars))

    # Add group-specific attributes
    ds.attrs.update(
        {
            "group_id": "auxiliary",
            "description": (
                "GNSS E/N/U displacement over the pair at coarse resolution: the"
                " UNR grid the calibration used, sampled with the same options"
                " (velocity x pair span; divide by pair_span_years for velocity)"
                if gnss_sampled
                else "Zero placeholder: this grid type has no GNSS sampling path"
            ),
            "units": "meters",
            "reference_frame": "ENU (East-North-Up)",
            "resolution": "coarse",
            "Conventions": "CF-1.8",
            "title": "Auxiliary 3D Displacement Model",
        }
    )
    if pair_span_years is not None:
        ds.attrs["pair_span_years"] = float(pair_span_years)

    return ds
