"""Build metadata group for calibration products."""

import numpy as np
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
    gnss_provenance_json: str | None = None,
    algorithm_schema_version: int | None = None,
    frame_parameters_version: str | None = None,
    components_applied_json: str | None = None,
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
    gnss_provenance_json : str or None, optional
        GNSS realisation record (Venti ``Provenance.as_dict()`` as JSON), written
        as ``gnss_provenance`` when given (v0.5 route); VLM verifies it
        (``venti.gnss.vlm``). Omitted otherwise, so the gamma product is
        unchanged.
    algorithm_schema_version : int or None, optional
        Venti algorithm-parameters schema the run used (v0.5 route).
    frame_parameters_version : str or None, optional
        Version of the frame-parameter table the run applied (v0.5 route,
        when the runconfig gives ``algorithm_parameters_overrides_json``).
    components_applied_json : str or None, optional
        Which ``cal_*`` components are in ``calibration`` (JSON, v0.5 route).

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

    if gnss_provenance_json is not None:
        data_vars["gnss_provenance"] = xr.DataArray(
            gnss_provenance_json,
            attrs={
                "description": (
                    "GNSS realisation the calibration was tied to (grid version and"
                    " type, reference frame, snapshot, lookup hash, buffer, defo DB"
                    " version, config and field hashes), JSON"
                ),
                "long_name": "GNSS Provenance",
                "dtype": "str",
            },
        )

    if algorithm_schema_version is not None:
        data_vars["algorithm_schema_version"] = xr.DataArray(
            np.int32(algorithm_schema_version),
            attrs={
                "description": "Venti algorithm-parameters schema version of the run",
                "long_name": "Algorithm Schema Version",
            },
        )
    if frame_parameters_version is not None:
        data_vars["frame_parameters_version"] = xr.DataArray(
            frame_parameters_version,
            attrs={
                "description": (
                    "Version of the frame-parameter table (per-frame algorithm"
                    " overrides) the run applied"
                ),
                "long_name": "Frame Parameters Version",
                "dtype": "str",
            },
        )
    if components_applied_json is not None:
        data_vars["components_applied"] = xr.DataArray(
            components_applied_json,
            attrs={
                "description": (
                    "Components summed into calibration (cal_gnss_surface,"
                    " cal_reference_offset, cal_tropo, cal_set, cal_unwrap_shift):"
                    " true when applied, JSON"
                ),
                "long_name": "Components Applied",
                "dtype": "str",
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
