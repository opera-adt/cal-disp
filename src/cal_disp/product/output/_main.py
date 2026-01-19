"""Build main group for calibration products."""

from datetime import datetime

import xarray as xr

from .._disp import DispProduct


def build_main_dataset(
    calibration: xr.DataArray,
    calibration_std: xr.DataArray | None,
    spatial_ref: xr.DataArray | None,
    disp_product: DispProduct,
    sensor: str,
    version: str,
    production_date: datetime,
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
    disp_product : DispProduct
        Original DISP product (for metadata).
    sensor : str
        Sensor type: "S1" or "NI".
    version : str
        Product version.
    production_date : datetime
        Date when product was generated.
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
        "product_type": f"OPERA_L4_CAL-DISP-{sensor}",
        "sensor": sensor,
        "frame_id": str(disp_product.frame_id),
        "mode": disp_product.mode,
        "polarization": disp_product.polarization,
        "reference_datetime": disp_product.reference_date.isoformat(),
        "secondary_datetime": disp_product.secondary_date.isoformat(),
        "production_datetime": production_date.isoformat(),
        "product_version": version,
        "group_id": "root",
        "description": (
            f"Calibration for {sensor} InSAR displacement (subtract from DISP)"
        ),
        "source_product": disp_product.filename,
        "usage": (
            "Subtract calibration layer from DISP displacement to obtain "
            "calibrated displacement"
        ),
        "Conventions": "CF-1.8",
        "title": f"OPERA L4 CAL-DISP-{sensor} Calibration Product",
        "institution": "NASA JPL",
        "contact": "operaops@jpl.nasa.gov",
    }

    ds.attrs.update(base_attrs)

    if metadata:
        ds.attrs.update(metadata)

    return ds
