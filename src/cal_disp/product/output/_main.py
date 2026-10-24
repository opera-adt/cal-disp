"""Build main group for calibration products."""

from __future__ import annotations

from datetime import datetime, timezone

import xarray as xr

from cal_disp._version import __version__

from ._utils import make_spatial_ref

# CF attributes of the projected coordinates, as in the DISP-S1 input
COORD_ATTRS: dict[str, dict[str, str]] = {
    "x": {
        "standard_name": "projection_x_coordinate",
        "long_name": "x coordinate of projection",
        "units": "m",
    },
    "y": {
        "standard_name": "projection_y_coordinate",
        "long_name": "y coordinate of projection",
        "units": "m",
    },
    "time": {
        "standard_name": "time",
        "long_name": "Time corresponding to beginning of secondary acquisition",
    },
}


def set_coord_attrs(ds: xr.Dataset) -> xr.Dataset:
    """Add the CF ``standard_name``/``long_name``/``units`` to ``x``/``y``/``time``.

    Attributes already present (e.g. copied from the DISP input) are kept.
    """
    for name, attrs in COORD_ATTRS.items():
        if name in ds.coords:
            for key, value in attrs.items():
                ds[name].attrs.setdefault(key, value)
    return ds


def build_main_dataset(
    calibration: xr.DataArray,
    calibration_std: xr.DataArray | None,
    spatial_ref: xr.DataArray | None,
    sensor: str,
    metadata: dict[str, str] | None,
    reference_date: datetime | None = None,
    secondary_date: datetime | None = None,
) -> xr.Dataset:
    """Build main dataset (root group) with calibration data.

    Parameters
    ----------
    calibration : xr.DataArray
        Calibration correction at full DISP resolution.
    calibration_std : xr.DataArray or None
        Calibration uncertainty at full resolution.
    spatial_ref : xr.DataArray or None
        Spatial reference data variable from input DISP product. Its CRS
        attributes are kept verbatim; ``GeoTransform`` is set for this grid.
    sensor : str
        Sensor type: "S1" or "NI".
    metadata : dict[str, str] or None
        Additional metadata.
    reference_date, secondary_date : datetime, optional
        Acquisition dates of the pair; give the ``temporal_resolution``
        attribute (the pair's temporal baseline).

    Returns
    -------
    xr.Dataset
        Main dataset with calibration data and attributes.

    """
    data_vars: dict[str, xr.DataArray] = {}

    # Calibration with description
    calibration = calibration.copy()
    calibration.attrs.update(
        {
            "description": "Calibration layer for DISP displacement",
            "long_name": "Calibration for DISP",
            "units": "meters",
            "grid_mapping": "spatial_ref",
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
            }
        )
        data_vars["calibration_std"] = calibration_std

    # Spatial reference: DISP CRS attributes + GeoTransform of this grid
    if spatial_ref is not None:
        data_vars["spatial_ref"] = make_spatial_ref(
            spatial_ref, calibration.x.values, calibration.y.values
        )

    ds = set_coord_attrs(xr.Dataset(data_vars))

    x = calibration.x.values
    spacing = float(abs(x[1] - x[0])) if x.size > 1 else float("nan")
    if reference_date is not None and secondary_date is not None:
        temporal_resolution = f"{(secondary_date - reference_date).days} days"
    else:
        temporal_resolution = "unknown"
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    base_attrs = {
        "Conventions": "CF-1.8",
        "title": f"OPERA L4 DISP-CAL-{sensor} Calibration Product",
        "institution": "NASA JPL",
        "contact": "operaops@jpl.nasa.gov",
        "source": "OPERA",
        "platform": sensor,
        "spatial_resolution": f"{spacing:g} meters",
        "temporal_resolution": temporal_resolution,
        "source_url": "https://www.jpl.nasa.gov/go/opera/products/disp-product-suite/",
        "references": "https://opera-adt.github.io/cal-disp/",
        "reference_document": "https://opera-adt.github.io/cal-disp/",
        "mission_name": "OPERA",
        "description": f"OPERA Calibration for {sensor} Surface Displacement product",
        "comment": (
            "Subtract calibration layer from DISP displacement to obtain "
            "calibrated displacement"
        ),
        # Software versions live in /metadata (cal_disp_software_version, ...)
        "history": f"{now}: created by cal_disp {__version__}",
    }

    ds.attrs.update(base_attrs)

    if metadata:
        ds.attrs.update(metadata)

    return ds
