"""Build auxiliary group for calibration products."""

import xarray as xr


def build_auxiliary_dataset(
    model_3d: dict[str, xr.DataArray] | None,
    model_3d_std: dict[str, xr.DataArray] | None,
    spatial_ref: xr.DataArray | None,
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
        Spatial reference data variable.

    Returns
    -------
    xr.Dataset
        Auxiliary dataset with 3D displacement model.

    """
    import rioxarray  # noqa: F401

    data_vars: dict[str, xr.DataArray] = {}

    # Component descriptions
    component_descriptions = {
        "north_south": "North-south displacement component (positive = north)",
        "east_west": "East-west displacement component (positive = east)",
        "up_down": "Up-down displacement component (positive = up)",
    }

    std_descriptions = {
        "north_south_std": "Uncertainty in north-south displacement",
        "east_west_std": "Uncertainty in east-west displacement",
        "up_down_std": "Uncertainty in up-down displacement",
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

    # Add spatial_ref to auxiliary group
    if spatial_ref is not None:
        data_vars["spatial_ref"] = spatial_ref

    ds = xr.Dataset(data_vars)

    # Write CRS if available
    if spatial_ref is not None:
        crs_wkt = spatial_ref.attrs.get("crs_wkt")
        if crs_wkt:
            ds = ds.rio.write_crs(crs_wkt)

    # Add group-specific attributes
    ds.attrs.update(
        {
            "group_id": "auxiliary",
            "description": (
                "3D displacement model at coarse resolution "
                "(e.g., from GNSS grid or deformation model)"
            ),
            "units": "meters",
            "reference_frame": "ENU (East-North-Up)",
            "resolution": "coarse",
        }
    )

    return ds
