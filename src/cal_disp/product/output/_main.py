"""Build main group for calibration products."""

import xarray as xr


def build_main_dataset(
    calibration: xr.DataArray,
    calibration_std: xr.DataArray | None,
    spatial_ref: xr.DataArray | None,
    sensor: str,
    metadata: dict[str, str] | None,
) -> xr.Dataset:
    """Build main dataset (root group) with calibration data.

    Parameters
    ----------
    calibration : xr.DataArray
        Calibration correction at full DISP resolution.
    calibration_std : xr.DataArray or None
        Calibration uncertainty at full resolution.
    spatial_ref : xr.DataArray or None
        Spatial reference data variable from input DISP product.
    sensor : str
        Sensor type: "S1" or "NI".
    metadata : dict[str, str] or None
        Additional metadata.

    Returns
    -------
    xr.Dataset
        Main dataset with calibration data and attributes.

    """
    import rioxarray  # noqa: F401

    data_vars: dict[str, xr.DataArray] = {}

    # Calibration with description
    calibration = calibration.copy()
    calibration.attrs.update(
        {
            "description": "Calibration layer for DISP displacement",
            "long_name": "Calibration for DISP",
            "units": "meters",
            "grid_mapping": "spatial_ref",
            "dtype": "float32",
            "coordinates": "time y x",
        }
    )
    data_vars["calibration"] = calibration

    # Calibration uncertainty
    if calibration_std is not None:
        calibration_std = calibration_std.copy()
        calibration_std.attrs.update(
            {
                "description": "Uncertainty in DISP calibration",
                "long_name": "DISP Calibration Uncertainty",
                "units": "meters",
                "grid_mapping": "spatial_ref",
                "dtype": "float32",
                "coordinates": "time y x",
            }
        )
        data_vars["calibration_std"] = calibration_std

    # Spatial reference
    if spatial_ref is not None:
        data_vars["spatial_ref"] = spatial_ref

    ds = xr.Dataset(data_vars)

    # Write CRS if available
    if spatial_ref is not None:
        crs_wkt = spatial_ref.attrs.get("crs_wkt")
        if crs_wkt:
            ds = ds.rio.write_crs(crs_wkt)

    # Add global attributes with type information
    base_attrs = {
        "Conventions": "CF-1.8",
        "title": f"OPERA L4 CAL-DISP-{sensor} Calibration Product",
        "institution": "NASA JPL",
        "contact": "operaops@jpl.nasa.gov",
        "source": "OPERA",
        "platform": sensor,
        "spatial_resolution": "30 meteres",
        "temporal_resolution": "12 days",
        "source_url": "https://www.jpl.nasa.gov/go/opera/products/disp-product-suite/",
        "references": "https://opera-adt.github.io/cal-disp/",
        "mision_name": "OPERA",
        "description": f"OPERA Calibration for {sensor} Surface Displacement product",
        "comment": (
            "Subtract calibration layer from DISP displacement to obtain "
            "calibrated displacement"
        ),
        "software": "cal_disp",
        "software_version": "0.1",
        "reference_document": "TBD",
        "history": "TBD",
    }

    ds.attrs.update(base_attrs)

    if metadata:
        ds.attrs.update(metadata)

    return ds
