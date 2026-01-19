"""Build metadata group for calibration products."""

import xarray as xr


def build_metadata_dataset(
    algorithm_parameters_yaml: str,
    platform_id: str,
    source_data_software_disp_version: str,
    cal_disp_software_version: str,
    venti_software_version: str,
    product_pixel_coordinate_convention: str = "center",
    ceos_atmospheric_phase_correction: str = "none",
    ceos_gridding_convention: str = "consistent",
    ceos_product_measurement_projection: str = "line_of_sight",
    ceos_ionospheric_phase_correction: str = "none",
    ceos_noise_removal: str = "N",
    algorithm_theoretical_basis_document_doi: str | None = None,
    product_landing_page_doi: str | None = None,
    product_specification_document_id: str | None = None,
    pge_runconfig: str | None = None,
) -> xr.Dataset:
    """Build metadata dataset with algorithm and processing information.

    Parameters
    ----------
    algorithm_parameters_yaml : str
        YAML string of algorithm parameters used.
    platform_id : str
        Sensor platform identification (e.g., "S1A", "S1B", "NISAR").
    source_data_software_disp_version : str
        DISP-S1 or DISP-NI processor version.
    cal_disp_software_version : str
        Version of cal_disp software.
    venti_software_version : str
        Version of venti software.
    product_pixel_coordinate_convention : str, optional
        Pixel coordinate convention: "center" or "corner". Default is "center".
    ceos_atmospheric_phase_correction : str, optional
        Atmospheric correction method. Default is "none".
    ceos_gridding_convention : str, optional
        Gridding/sampling frame consistency. Default is "consistent".
    ceos_product_measurement_projection : str, optional
        Displacement projection: "line_of_sight", "horizontal", or "vertical".
        Default is "line_of_sight".
    ceos_ionospheric_phase_correction : str, optional
        Ionospheric correction method. Default is "none".
    ceos_noise_removal : str, optional
        Noise removal flag: "Y" or "N". Default is "N".
    algorithm_theoretical_basis_document_doi : str or None, optional
        DOI for algorithm theoretical basis document.
    product_landing_page_doi : str or None, optional
        DOI for product landing page.
    product_specification_document_id : str or None, optional
        Product specification document identifier.
    pge_runconfig : str or None, optional
        Full PGE runconfig YAML file content.

    Returns
    -------
    xr.Dataset
        Metadata dataset with scalar variables.

    """
    data_vars: dict[str, xr.DataArray] = {}

    # Algorithm parameters
    data_vars["algorithm_parameters_yaml"] = xr.DataArray(
        algorithm_parameters_yaml,
        attrs={
            "description": "Algorithm parameters in YAML format",
            "long_name": "Algorithm Parameters",
            "dtype": "str",
            "format": "YAML",
        },
    )

    # Platform ID
    data_vars["platform_id"] = xr.DataArray(
        platform_id,
        attrs={
            "description": "Sensor platform identification string",
            "long_name": "Platform ID",
            "examples": "S1A, S1B, NISAR",
            "dtype": "str",
        },
    )

    # Source data software version
    data_vars["source_data_software_disp_version"] = xr.DataArray(
        source_data_software_disp_version,
        attrs={
            "description": "DISP-S1 or DISP-NI processor version used for source data",
            "long_name": "Source DISP Software Version",
            "dtype": "str",
        },
    )

    # cal_disp software version
    data_vars["cal_disp_software_version"] = xr.DataArray(
        cal_disp_software_version,
        attrs={
            "description": "Version of cal_disp software used to generate product",
            "long_name": "cal_disp Software Version",
            "dtype": "str",
        },
    )

    # venti software version
    data_vars["venti_software_version"] = xr.DataArray(
        venti_software_version,
        attrs={
            "description": "Version of venti software used to generate product",
            "long_name": "Venti Software Version",
            "dtype": "str",
        },
    )

    # Pixel coordinate convention
    data_vars["product_pixel_coordinate_convention"] = xr.DataArray(
        product_pixel_coordinate_convention,
        attrs={
            "description": (
                "x/y coordinate convention referring to pixel center or corner"
            ),
            "long_name": "Pixel Coordinate Convention",
            "valid_values": "center, corner",
            "dtype": "str",
        },
    )

    # CEOS metadata fields
    data_vars["ceos_atmospheric_phase_correction"] = xr.DataArray(
        ceos_atmospheric_phase_correction,
        attrs={
            "description": "Method used to correct for atmospheric phase noise",
            "long_name": "Atmospheric Phase Correction",
            "dtype": "str",
        },
    )

    data_vars["ceos_gridding_convention"] = xr.DataArray(
        ceos_gridding_convention,
        attrs={
            "description": (
                "Whether a consistent gridding/sampling frame is used for "
                "ascending/descending frames"
            ),
            "long_name": "Gridding Convention",
            "dtype": "str",
        },
    )

    data_vars["ceos_product_measurement_projection"] = xr.DataArray(
        ceos_product_measurement_projection,
        attrs={
            "description": "Projection of the displacement image",
            "long_name": "Product Measurement Projection",
            "valid_values": "line_of_sight, horizontal, vertical",
            "dtype": "str",
        },
    )

    data_vars["ceos_ionospheric_phase_correction"] = xr.DataArray(
        ceos_ionospheric_phase_correction,
        attrs={
            "description": "Method used to correct for ionospheric phase noise",
            "long_name": "Ionospheric Phase Correction",
            "dtype": "str",
        },
    )

    data_vars["ceos_noise_removal"] = xr.DataArray(
        ceos_noise_removal,
        attrs={
            "description": (
                "Flag if noise removal has been applied. Should include noise removal "
                "algorithm and reference (URL or DOI) if Y"
            ),
            "long_name": "Noise Removal Flag",
            "valid_values": "Y, N",
            "dtype": "str",
        },
    )

    # Optional DOIs and document IDs
    if algorithm_theoretical_basis_document_doi is not None:
        data_vars["algorithm_theoretical_basis_document_doi"] = xr.DataArray(
            algorithm_theoretical_basis_document_doi,
            attrs={
                "description": "DOI for algorithm theoretical basis document",
                "long_name": "Algorithm Theoretical Basis Document DOI",
                "dtype": "str",
                "format": "DOI",
            },
        )

    if product_landing_page_doi is not None:
        data_vars["product_landing_page_doi"] = xr.DataArray(
            product_landing_page_doi,
            attrs={
                "description": "DOI for product landing page",
                "long_name": "Product Landing Page DOI",
                "dtype": "str",
                "format": "DOI",
            },
        )

    if product_specification_document_id is not None:
        data_vars["product_specification_document_id"] = xr.DataArray(
            product_specification_document_id,
            attrs={
                "description": "Document identifier for Product Specification",
                "long_name": "Product Specification Document ID",
                "dtype": "str",
            },
        )

    # Optional PGE runconfig
    if pge_runconfig is not None:
        data_vars["pge_runconfig"] = xr.DataArray(
            pge_runconfig,
            attrs={
                "description": "Full PGE runconfig YAML file used to generate product",
                "long_name": "PGE Run Configuration",
                "dtype": "str",
                "format": "YAML",
            },
        )

    ds = xr.Dataset(data_vars)

    # Add group attributes
    ds.attrs.update(
        {
            "group_id": "metadata",
            "description": (
                "Contains algorithm parameters, processing metadata, and CEOS-compliant"
                " metadata fields."
            ),
        }
    )

    return ds
