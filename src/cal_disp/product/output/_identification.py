"""Build identification group for calibration products."""

from datetime import datetime

import xarray as xr


def build_identification_dataset(
    frame_id: int,
    reference_datetime: datetime,
    secondary_datetime: datetime,
    product_version: str,
    calibration_reference_name: str,
    calibration_reference_version: str,
    calibration_reference_type: str,
    calibration_reference_reference_frame: str,
    absolute_orbit_number: int,
    acquisition_mode: str,
    bounding_polygon: str,
    instrument_name: str,
    look_direction: str,
    nodata_pixel_count: int,
    orbit_pass_direction: str,
    processing_facility: str,
    processing_start_datetime: datetime,
    product_bounding_box: str,
    product_data_access: str,
    product_data_polarization: str,
    product_sample_spacing: str,
    radar_band: str,
    source_data_access: str,
    source_data_dem_name: str,
    source_data_file_list: list[str],
    source_calibration_file_list: list[str],
    source_data_imaging_geometry: str,
    source_data_satellite_names: list[str],
    source_data_x_spacing: float,
    source_data_y_spacing: float,
    static_layers_data_access: str,
    track_number: int,
    ceos_number_of_input_granules: int,
    ceos_analysis_ready_data_document_identifier: str | None = None,
    ceos_analysis_ready_data_product_type: str | None = None,
) -> xr.Dataset:
    """Build identification dataset with product metadata.

    Parameters
    ----------
    frame_id : int
        OPERA frame identifier.
    reference_datetime : datetime
        Reference (earlier) acquisition datetime.
    secondary_datetime : datetime
        Secondary (later) acquisition datetime.
    product_version : str
        Product specification version.
    calibration_reference_name : str
        Name of calibration reference data source.
    calibration_reference_version : str
        Version of calibration reference data.
    calibration_reference_type : str
        Type of calibration reference data (e.g., "constant", "variable").
    calibration_reference_reference_frame : str
        Reference frame of calibration data (e.g., "IGS20").
    absolute_orbit_number : int
        Absolute orbit number.
    acquisition_mode : str
        Radar acquisition mode for input products.
    bounding_polygon : str
        WKT representation of bounding polygon.
    instrument_name : str
        Instrument name (e.g., "C-SAR").
    look_direction : str
        Look direction: "left" or "right".
    nodata_pixel_count : int
        Number of nodata pixels.
    orbit_pass_direction : str
        Orbit direction: "ascending" or "descending".
    processing_facility : str
        Product processing facility.
    processing_start_datetime : datetime
        UTC datetime of processing start.
    product_bounding_box : str
        Opposite corners in UTM coordinates as (west, south, east, north).
    product_data_access : str
        URL or DOI for product data access.
    product_data_polarization : str
        Radar polarization of displacement products.
    product_sample_spacing : str
        Spacing between adjacent X/Y samples in UTM coordinates.
    radar_band : str
        Acquired radar frequency band (e.g., "C", "L").
    source_data_access : str
        URL or DOI for source data access.
    source_data_dem_name : str
        Name of DEM used during input data processing.
    source_data_file_list : list[str]
        List of input DISP granules used to create calibration layer.
    source_calibration_file_list : list[str]
        List of input UNR latlon and grid points tenv8 files.
    source_data_imaging_geometry : str
        Imaging geometry of DISP.
    source_data_satellite_names : list[str]
        Names of satellites included in input granules.
    source_data_x_spacing : float
        Pixel spacing of source geocoded SLC data in x-direction.
    source_data_y_spacing : float
        Pixel spacing of source geocoded SLC data in y-direction.
    static_layers_data_access : str
        URL of static layers product associated with this calibration product.
    track_number : int
        Track number.
    ceos_number_of_input_granules : int
        Number of input data granules used during processing.
    ceos_analysis_ready_data_document_identifier : str or None, optional
        CEOS CARD document identifier.
    ceos_analysis_ready_data_product_type : str or None, optional
        CEOS CARD product type name.

    Returns
    -------
    xr.Dataset
        Identification dataset with scalar variables.

    """
    data_vars: dict[str, xr.DataArray] = {}

    # Product identification
    data_vars["frame_id"] = xr.DataArray(
        frame_id,
        attrs={
            "description": "OPERA frame identifier",
            "long_name": "Frame ID",
            "dtype": "int32",
        },
    )

    data_vars["product_version"] = xr.DataArray(
        product_version,
        attrs={
            "description": "Product specification version",
            "long_name": "Product Version",
            "dtype": "str",
        },
    )

    data_vars["track_number"] = xr.DataArray(
        track_number,
        attrs={
            "description": "Track number",
            "long_name": "Track Number",
            "dtype": "int32",
        },
    )

    data_vars["absolute_orbit_number"] = xr.DataArray(
        absolute_orbit_number,
        attrs={
            "description": "Absolute orbit number",
            "long_name": "Absolute Orbit Number",
            "dtype": "int32",
        },
    )

    data_vars["orbit_pass_direction"] = xr.DataArray(
        orbit_pass_direction,
        attrs={
            "description": "Orbit direction: ascending or descending",
            "long_name": "Orbit Pass Direction",
            "valid_values": "ascending, descending",
            "dtype": "str",
        },
    )

    # Temporal information
    data_vars["reference_datetime"] = xr.DataArray(
        reference_datetime.isoformat(),
        attrs={
            "description": (
                "Reference (earlier) acquisition datetime in ISO 8601 format"
            ),
            "long_name": "Reference Datetime",
            "dtype": "str",
            "format": "ISO 8601",
        },
    )

    data_vars["secondary_datetime"] = xr.DataArray(
        secondary_datetime.isoformat(),
        attrs={
            "description": "Secondary (later) acquisition datetime in ISO 8601 format",
            "long_name": "Secondary Datetime",
            "dtype": "str",
            "format": "ISO 8601",
        },
    )

    data_vars["processing_start_datetime"] = xr.DataArray(
        processing_start_datetime.isoformat(),
        attrs={
            "description": "UTC datetime of the start of processing for this product",
            "long_name": "Processing Start Datetime",
            "dtype": "str",
            "format": "ISO 8601",
        },
    )

    # Spatial information
    data_vars["bounding_polygon"] = xr.DataArray(
        bounding_polygon,
        attrs={
            "description": "WKT representation of bounding polygon of the image",
            "long_name": "Bounding Polygon",
            "dtype": "str",
            "format": "WKT",
        },
    )

    data_vars["product_bounding_box"] = xr.DataArray(
        product_bounding_box,
        attrs={
            "description": (
                "Opposite corners in UTM coordinates as (west, south, east, north)"
            ),
            "long_name": "Product Bounding Box",
            "dtype": "str",
        },
    )

    data_vars["product_sample_spacing"] = xr.DataArray(
        product_sample_spacing,
        attrs={
            "description": "Spacing between adjacent X/Y samples in UTM coordinates",
            "long_name": "Product Sample Spacing",
            "dtype": "str",
        },
    )

    # Calibration reference information
    data_vars["calibration_reference_name"] = xr.DataArray(
        calibration_reference_name,
        attrs={
            "description": "Name of calibration reference data source",
            "long_name": "Calibration Reference Name",
            "dtype": "str",
        },
    )

    data_vars["calibration_reference_version"] = xr.DataArray(
        calibration_reference_version,
        attrs={
            "description": "Version of calibration reference data",
            "long_name": "Calibration Reference Version",
            "dtype": "str",
        },
    )

    data_vars["calibration_reference_type"] = xr.DataArray(
        calibration_reference_type,
        attrs={
            "description": "Type of calibration reference data (constant or variable)",
            "long_name": "Calibration Reference Type",
            "valid_values": "constant, variable",
            "dtype": "str",
        },
    )

    data_vars["calibration_reference_reference_frame"] = xr.DataArray(
        calibration_reference_reference_frame,
        attrs={
            "description": "Reference frame of calibration data (e.g., IGS20)",
            "long_name": "Calibration Reference Frame",
            "dtype": "str",
        },
    )

    data_vars["source_calibration_file_list"] = xr.DataArray(
        ", ".join(source_calibration_file_list),
        attrs={
            "description": (
                "Comma-separated list of input UNR latlon and grid points tenv8 files"
            ),
            "long_name": "Source Calibration File List",
            "dtype": "str",
        },
    )

    # Source data information
    data_vars["source_data_file_list"] = xr.DataArray(
        ", ".join(source_data_file_list),
        attrs={
            "description": (
                "Comma-separated list of input DISP granules used to create calibration"
                " layer"
            ),
            "long_name": "Source Data File List",
            "dtype": "str",
        },
    )

    data_vars["source_data_access"] = xr.DataArray(
        source_data_access,
        attrs={
            "description": "URL or DOI for source data retrieval",
            "long_name": "Source Data Access",
            "dtype": "str",
        },
    )

    data_vars["source_data_dem_name"] = xr.DataArray(
        source_data_dem_name,
        attrs={
            "description": "Name of DEM used during input data processing",
            "long_name": "Source Data DEM Name",
            "dtype": "str",
        },
    )

    data_vars["source_data_satellite_names"] = xr.DataArray(
        ", ".join(source_data_satellite_names),
        attrs={
            "description": (
                "Comma-separated names of satellites included in input granules"
            ),
            "long_name": "Source Data Satellite Names",
            "dtype": "str",
        },
    )

    data_vars["source_data_imaging_geometry"] = xr.DataArray(
        source_data_imaging_geometry,
        attrs={
            "description": "Imaging geometry of DISP",
            "long_name": "Source Data Imaging Geometry",
            "dtype": "str",
        },
    )

    data_vars["source_data_x_spacing"] = xr.DataArray(
        source_data_x_spacing,
        attrs={
            "description": "Pixel spacing of source geocoded SLC data in x-direction",
            "long_name": "Source Data X Spacing",
            "units": "meters",
            "dtype": "float64",
        },
    )

    data_vars["source_data_y_spacing"] = xr.DataArray(
        source_data_y_spacing,
        attrs={
            "description": "Pixel spacing of source geocoded SLC data in y-direction",
            "long_name": "Source Data Y Spacing",
            "units": "meters",
            "dtype": "float64",
        },
    )

    # Radar configuration
    data_vars["acquisition_mode"] = xr.DataArray(
        acquisition_mode,
        attrs={
            "description": "Radar acquisition mode for input products",
            "long_name": "Acquisition Mode",
            "dtype": "str",
        },
    )

    data_vars["instrument_name"] = xr.DataArray(
        instrument_name,
        attrs={
            "description": "Instrument name",
            "long_name": "Instrument Name",
            "dtype": "str",
        },
    )

    data_vars["look_direction"] = xr.DataArray(
        look_direction,
        attrs={
            "description": "Look direction: left or right",
            "long_name": "Look Direction",
            "valid_values": "left, right",
            "dtype": "str",
        },
    )

    data_vars["product_data_polarization"] = xr.DataArray(
        product_data_polarization,
        attrs={
            "description": "Radar polarization of displacement products",
            "long_name": "Product Data Polarization",
            "dtype": "str",
        },
    )

    data_vars["radar_band"] = xr.DataArray(
        radar_band,
        attrs={
            "description": "Acquired radar frequency band",
            "long_name": "Radar Band",
            "dtype": "str",
        },
    )

    # Processing information
    data_vars["processing_facility"] = xr.DataArray(
        processing_facility,
        attrs={
            "description": "Product processing facility",
            "long_name": "Processing Facility",
            "dtype": "str",
        },
    )

    data_vars["ceos_number_of_input_granules"] = xr.DataArray(
        ceos_number_of_input_granules,
        attrs={
            "description": "Number of input data granules used during processing",
            "long_name": "Number of Input Granules",
            "dtype": "int32",
        },
    )

    data_vars["nodata_pixel_count"] = xr.DataArray(
        nodata_pixel_count,
        attrs={
            "description": "Number of nodata pixels",
            "long_name": "NoData Pixel Count",
            "dtype": "int64",
        },
    )

    # Data access
    data_vars["product_data_access"] = xr.DataArray(
        product_data_access,
        attrs={
            "description": "URL or DOI for product data access",
            "long_name": "Product Data Access",
            "dtype": "str",
        },
    )

    data_vars["static_layers_data_access"] = xr.DataArray(
        static_layers_data_access,
        attrs={
            "description": (
                "URL of static layers product associated with this calibration product"
            ),
            "long_name": "Static Layers Data Access",
            "dtype": "str",
        },
    )

    # Optional CEOS fields
    if ceos_analysis_ready_data_document_identifier is not None:
        data_vars["ceos_analysis_ready_data_document_identifier"] = xr.DataArray(
            ceos_analysis_ready_data_document_identifier,
            attrs={
                "description": "CEOS Analysis Ready Data (CARD) document identifier",
                "long_name": "CEOS CARD Document Identifier",
                "dtype": "str",
            },
        )

    if ceos_analysis_ready_data_product_type is not None:
        data_vars["ceos_analysis_ready_data_product_type"] = xr.DataArray(
            ceos_analysis_ready_data_product_type,
            attrs={
                "description": "CEOS Analysis Ready Data (CARD) product type name",
                "long_name": "CEOS CARD Product Type",
                "dtype": "str",
            },
        )

    ds = xr.Dataset(data_vars)

    # Add group attributes
    ds.attrs.update(
        {
            "group_id": "identification",
            "description": (
                "Contains comprehensive product identification metadata including "
                "temporal, spatial, calibration source, and processing information."
            ),
        }
    )

    return ds
