from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
import xarray as xr
from rasterio.crs import CRS
from rasterio.warp import transform_bounds

from .._disp import DispProduct
from ._auxiliary import build_auxiliary_dataset
from ._identification import build_identification_dataset
from ._main import build_main_dataset
from ._metadata import build_metadata_dataset
from ._utils import (
    build_filename,
    compute_transform_from_coords,
    get_crs,
    get_transform,
)


@dataclass
class CalProduct:
    """Calibrated OPERA DISP displacement product.

    Represents a calibration correction product that should be subtracted
    from OPERA DISP products. Main group contains calibration at full
    resolution. Optional auxiliary group contains 3D displacement
    components at coarser resolution. Identification group contains
    product metadata and provenance.

    Groups
    ------
    Main group:
        - calibration: Correction to subtract from DISP (full resolution)
        - calibration_std: Calibration uncertainty (full resolution)

    Identification group:
        - frame_id: OPERA frame identifier
        - product_version: Product version string
        - reference_datetime: Reference acquisition datetime
        - secondary_datetime: Secondary acquisition datetime
        - source_calibration_reference_name: Calibration source name
        - source_calibration_reference_version: Calibration source version
        - source_calibration_reference_type: Calibration source type
        - source_calibration_reference_reference_frame: Reference frame
        - pge_runconfig: PGE run configuration (optional)
        - software_version: Software version (optional)

    Metadata group (optional):
        - algorithm_parameters_yaml, platform_id
        - Software versions (cal_disp, venti, DISP)
        - CEOS compliance fields
        - Product specification information

    Auxiliary group (optional):
        - north_south: North-south displacement component (coarse resolution)
        - east_west: East-west displacement component (coarse resolution)
        - up_down: Up-down displacement component (coarse resolution)
        - north_south_std: Uncertainty in north-south
        - east_west_std: Uncertainty in east-west
        - up_down_std: Uncertainty in up-down

    Parameters
    ----------
    path : Path
        Path to the calibration product NetCDF file.
    frame_id : int
        OPERA frame identifier.
    primary_date : datetime
        Earlier acquisition date (reference).
    secondary_date : datetime
        Later acquisition date.
    polarization : str
        Radar polarization (e.g., "VV", "VH").
    sensor : str
        Sensor type: "S1" (Sentinel-1) or "NI" (NISAR).
    version : str
        Product version string.
    production_date : datetime
        Date when product was generated.
    mode : str
        Acquisition mode (e.g., "IW" for S1, "LSAR" for NI).

    Examples
    --------
    >>> # Create product with both groups
    >>> cal = CalProduct.create(
    ...     calibration=cal_correction,
    ...     disp_product=disp,
    ...     output_dir="output/",
    ... )

    >>> # Access main calibration
    >>> ds_main = cal.open_dataset()
    >>> calibration = ds_main["calibration"]

    >>> # Access auxiliary 3D model (coarse resolution)
    >>> ds_model = cal.open_auxiliary()
    >>> model_up = ds_model["up_down"]

    """

    path: Path
    frame_id: int
    reference_date: datetime
    secondary_date: datetime
    polarization: str
    sensor: str
    version: str
    production_date: datetime
    mode: str = "IW"

    # Filename pattern supporting both S1 and NI sensors
    _PATTERN = re.compile(
        r"OPERA_L4_CAL-DISP-(?P<sensor>S1|NI)_"
        r"(?P<mode>\w+)_"
        r"F(?P<frame_id>\d+)_"
        r"(?P<pol>\w+)_"
        r"(?P<reference>\d{8}T\d{6}Z)_"
        r"(?P<secondary>\d{8}T\d{6}Z)_"
        r"v(?P<version>[\d.]+)_"
        r"(?P<production>\d{8}T\d{6}Z)"
        r"\.nc$"
    )

    # Main group layers (full resolution)
    CAL_LAYERS = [
        "calibration",
        "calibration_std",
    ]

    # model_3d group layers (coarse resolution)
    MODEL_3D_LAYERS = [
        "north_south",
        "east_west",
        "up_down",
        "north_south_std",
        "east_west_std",
        "up_down_std",
    ]

    def __post_init__(self) -> None:
        """Validate product after construction."""
        self.path = Path(self.path)

        if self.frame_id <= 0:
            raise ValueError(f"frame_id must be positive, got {self.frame_id}")

        if self.secondary_date <= self.reference_date:
            raise ValueError(
                f"Secondary date ({self.secondary_date}) must be after "
                f"Reference date ({self.reference_date})"
            )

        if self.polarization not in {"VV", "VH", "HH", "HV"}:
            raise ValueError(f"Invalid polarization: {self.polarization}")

        if self.sensor not in {"S1", "NI"}:
            raise ValueError(f"Invalid sensor: {self.sensor}. Must be 'S1' or 'NI'")

    @classmethod
    def from_path(cls, path: Path | str) -> "CalProduct":
        """Parse product metadata from filename.

        Parameters
        ----------
        path : Path or str
            Path to calibration product NetCDF file.

        Returns
        -------
        CalProduct
            Parsed calibration product instance.

        Raises
        ------
        ValueError
            If filename doesn't match OPERA CAL-DISP pattern.

        Examples
        --------
        >>> cal = CalProduct.from_path(
        ...     "OPERA_L4_CAL-DISP-S1_IW_F08882_VV_20220111T002651Z_"
        ...     "20220722T002657Z_v1.0_20251227T123456Z.nc"
        ... )
        >>> cal.sensor
        'S1'

        """
        path = Path(path)
        match = cls._PATTERN.match(path.name)

        if not match:
            raise ValueError(
                f"Filename does not match OPERA CAL-DISP pattern: {path.name}"
            )

        return cls(
            path=path,
            frame_id=int(match.group("frame_id")),
            reference_date=datetime.strptime(
                match.group("reference"), "%Y%m%dT%H%M%SZ"
            ),
            secondary_date=datetime.strptime(
                match.group("secondary"), "%Y%m%dT%H%M%SZ"
            ),
            polarization=match.group("pol"),
            sensor=match.group("sensor"),
            version=match.group("version"),
            production_date=datetime.strptime(
                match.group("production"), "%Y%m%dT%H%M%SZ"
            ),
            mode=match.group("mode"),
        )

    @classmethod
    def create(
        cls,
        calibration: xr.DataArray,
        disp_product: DispProduct,
        output_dir: Path | str,
        sensor: str = "S1",
        calibration_std: xr.DataArray | None = None,
        spatial_ref: xr.DataArray | None = None,
        global_metadata: dict[str, str] | None = None,
        version: str = "1.0",
    ) -> CalProduct:
        """Create calibration product with main group only.

        Parameters
        ----------
        calibration : xr.DataArray
            Calibration correction with dimensions (time, y, x).
        disp_product : DispProduct
            Original DISP product (for metadata).
        output_dir : Path or str
            Output directory for NetCDF file.
        sensor : str, optional
            Sensor type: "S1" or "NI". Default is "S1".
        calibration_std : xr.DataArray or None, optional
            Calibration uncertainty. Default is None.
        spatial_ref : xr.DataArray or None, optional
            Spatial reference from input DISP product. Default is None.
        global_metadata : dict[str, str] or None, optional
            Additional metadata. Default is None.
        version : str, optional
            Product version. Default is "1.0".

        Returns
        -------
        CalProduct
            Created calibration product with main group only.

        Raises
        ------
        ValueError
            If sensor is not 'S1' or 'NI' or calibration missing time dimension.

        Examples
        --------
        >>> cal = CalProduct.create(
        ...     calibration=cal_correction,
        ...     disp_product=disp,
        ...     output_dir="output/",
        ...     calibration_std=cal_std,
        ...     spatial_ref=ds["spatial_ref"],
        ... )
        >>> cal.add_identification(unr_grid_version="2024.1")
        >>> cal.add_metadata(platform_id="S1A", ...)
        >>> cal.add_auxiliary(model_3d=model_dict)

        """
        # Validate inputs
        if sensor not in {"S1", "NI"}:
            raise ValueError(f"Invalid sensor: {sensor}. Must be 'S1' or 'NI'")

        if "time" not in calibration.dims:
            raise ValueError("calibration must have 'time' dimension")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate filename
        production_date = datetime.utcnow()
        filename = build_filename(
            disp_product,
            sensor,
            version,
            production_date,
        )
        output_file = output_dir / filename

        # Build and save main dataset
        ds_main = build_main_dataset(
            calibration=calibration,
            calibration_std=calibration_std,
            spatial_ref=spatial_ref,
            disp_product=disp_product,
            sensor=sensor,
            version=version,
            production_date=production_date,
            metadata=global_metadata,
        )
        ds_main.to_netcdf(output_file, engine="h5netcdf")

        return cls(
            path=output_file,
            frame_id=disp_product.frame_id,
            reference_date=disp_product.reference_date,
            secondary_date=disp_product.secondary_date,
            polarization=disp_product.polarization,
            sensor=sensor,
            version=version,
            production_date=production_date,
            mode=disp_product.mode,
        )

    def add_identification(
        self,
        calibration_reference_name: str,
        calibration_reference_version: str,
        calibration_reference_type: str,
        calibration_reference_reference_frame: str,
        source_data_file_list: list[str],
        source_calibration_file_list: list[str],
        source_data_access: str,
        source_data_dem_name: str,
        source_data_satellite_names: list[str],
        source_data_imaging_geometry: str,
        source_data_x_spacing: float,
        source_data_y_spacing: float,
        static_layers_data_access: str,
        absolute_orbit_number: int,
        track_number: int,
        instrument_name: str,
        look_direction: str,
        radar_band: str,
        orbit_pass_direction: str,
        bounding_polygon: str,
        product_bounding_box: str,
        product_sample_spacing: str,
        product_data_access: str,
        processing_facility: str,
        nodata_pixel_count: int,
        ceos_number_of_input_granules: int,
        processing_start_datetime: datetime | None = None,
        ceos_analysis_ready_data_document_identifier: str | None = None,
        ceos_analysis_ready_data_product_type: str | None = None,
    ) -> None:
        """Add identification group to existing product.

        Parameters
        ----------
        calibration_reference_name : str
            Name of calibration reference data source (e.g., "UNR gridded data").
        calibration_reference_version : str
            Version of calibration reference data.
        calibration_reference_type : str
            Type of calibration reference: "constant" or "variable".
        calibration_reference_reference_frame : str
            Reference frame (e.g., "IGS20").
        source_data_file_list : list[str]
            List of input DISP granules used.
        source_calibration_file_list : list[str]
            List of UNR calibration files used.
        source_data_access : str
            URL or DOI for source DISP data.
        source_data_dem_name : str
            DEM name used in processing.
        source_data_satellite_names : list[str]
            Satellite names (e.g., ["Sentinel-1A", "Sentinel-1B"]).
        source_data_imaging_geometry : str
            DISP imaging geometry.
        source_data_x_spacing : float
            Source data pixel spacing in x-direction (meters).
        source_data_y_spacing : float
            Source data pixel spacing in y-direction (meters).
        static_layers_data_access : str
            URL for static layers product.
        absolute_orbit_number : int
            Absolute orbit number.
        track_number : int
            Track number.
        instrument_name : str
            Instrument name (e.g., "C-SAR").
        look_direction : str
            Look direction: "left" or "right".
        radar_band : str
            Radar frequency band (e.g., "C", "L").
        orbit_pass_direction : str
            Orbit direction: "ascending" or "descending".
        bounding_polygon : str
            WKT bounding polygon.
        product_bounding_box : str
            Bounding box as "(west, south, east, north)".
        product_sample_spacing : str
            Product sample spacing in UTM coordinates.
        product_data_access : str
            URL or DOI for product access.
        processing_facility : str
            Processing facility name.
        nodata_pixel_count : int
            Number of nodata pixels.
        ceos_number_of_input_granules : int
            Number of input granules used.
        processing_start_datetime : datetime or None, optional
            Processing start time. Defaults to current UTC time.
        ceos_analysis_ready_data_document_identifier : str or None, optional
            CEOS CARD document identifier.
        ceos_analysis_ready_data_product_type : str or None, optional
            CEOS CARD product type.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.

        Examples
        --------
        >>> cal.add_identification(
        ...     calibration_reference_name="UNR gridded data",
        ...     calibration_reference_version="2024.1",
        ...     calibration_reference_type="constant",
        ...     calibration_reference_reference_frame="IGS20",
        ...     source_data_file_list=["DISP_file1.nc", "DISP_file2.nc"],
        ...     source_calibration_file_list=["unr_latlon.txt", "station1.tenv8"],
        ...     source_data_access="https://doi.org/...",
        ...     source_data_dem_name="Copernicus DEM GLO-30",
        ...     source_data_satellite_names=["Sentinel-1A"],
        ...     source_data_imaging_geometry="right_looking",
        ...     source_data_x_spacing=30.0,
        ...     source_data_y_spacing=30.0,
        ...     static_layers_data_access="https://...",
        ...     absolute_orbit_number=12345,
        ...     track_number=64,
        ...     instrument_name="C-SAR",
        ...     look_direction="right",
        ...     radar_band="C",
        ...     orbit_pass_direction="ascending",
        ...     bounding_polygon="POLYGON((...)",
        ...     product_bounding_box="(500000, 3500000, 600000, 3600000)",
        ...     product_sample_spacing="30m",
        ...     product_data_access="https://...",
        ...     processing_facility="JPL",
        ...     nodata_pixel_count=1234,
        ...     ceos_number_of_input_granules=2,
        ... )

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        # Default processing start time to now if not provided
        if processing_start_datetime is None:
            processing_start_datetime = datetime.utcnow()

        ds_id = build_identification_dataset(
            frame_id=self.frame_id,
            reference_datetime=self.reference_date,
            secondary_datetime=self.secondary_date,
            product_version=self.version,
            calibration_reference_name=calibration_reference_name,
            calibration_reference_version=calibration_reference_version,
            calibration_reference_type=calibration_reference_type,
            calibration_reference_reference_frame=calibration_reference_reference_frame,
            absolute_orbit_number=absolute_orbit_number,
            acquisition_mode=self.mode,  # From CalProduct
            bounding_polygon=bounding_polygon,
            instrument_name=instrument_name,
            look_direction=look_direction,
            nodata_pixel_count=nodata_pixel_count,
            orbit_pass_direction=orbit_pass_direction,
            processing_facility=processing_facility,
            processing_start_datetime=processing_start_datetime,
            product_bounding_box=product_bounding_box,
            product_data_access=product_data_access,
            product_data_polarization=self.polarization,  # From CalProduct
            product_sample_spacing=product_sample_spacing,
            radar_band=radar_band,
            source_data_access=source_data_access,
            source_data_dem_name=source_data_dem_name,
            source_data_file_list=source_data_file_list,
            source_calibration_file_list=source_calibration_file_list,
            source_data_imaging_geometry=source_data_imaging_geometry,
            source_data_satellite_names=source_data_satellite_names,
            source_data_x_spacing=source_data_x_spacing,
            source_data_y_spacing=source_data_y_spacing,
            static_layers_data_access=static_layers_data_access,
            track_number=track_number,
            ceos_number_of_input_granules=ceos_number_of_input_granules,
            ceos_analysis_ready_data_document_identifier=ceos_analysis_ready_data_document_identifier,
            ceos_analysis_ready_data_product_type=ceos_analysis_ready_data_product_type,
        )
        ds_id.to_netcdf(self.path, mode="a", group="identification", engine="h5netcdf")

    def add_metadata(
        self,
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
    ) -> None:
        """Add metadata group to existing product.

        Parameters
        ----------
        algorithm_parameters_yaml : str
            YAML string of algorithm parameters.
        platform_id : str
            Sensor platform ID (e.g., "S1A", "S1B").
        source_data_software_disp_version : str
            DISP-S1/DISP-NI processor version.
        cal_disp_software_version : str
            cal_disp software version.
        venti_software_version : str
            venti software version.
        product_pixel_coordinate_convention : str, optional
            Pixel coordinate convention. Default is "center".
        ceos_atmospheric_phase_correction : str, optional
            Atmospheric correction method. Default is "none".
        ceos_gridding_convention : str, optional
            Gridding convention. Default is "consistent".
        ceos_product_measurement_projection : str, optional
            Measurement projection. Default is "line_of_sight".
        ceos_ionospheric_phase_correction : str, optional
            Ionospheric correction method. Default is "none".
        ceos_noise_removal : str, optional
            Noise removal flag. Default is "N".
        algorithm_theoretical_basis_document_doi : str or None, optional
            Algorithm ATBD DOI.
        product_landing_page_doi : str or None, optional
            Product landing page DOI.
        product_specification_document_id : str or None, optional
            Product spec document ID.
        pge_runconfig : str or None, optional
            PGE runconfig YAML content.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.

        Examples
        --------
        >>> cal.add_metadata(
        ...     algorithm_parameters_yaml=yaml_str,
        ...     platform_id="S1A",
        ...     source_data_software_disp_version="1.0.0",
        ...     cal_disp_software_version="0.1.0",
        ...     venti_software_version="0.2.0",
        ... )

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        ds_meta = build_metadata_dataset(
            algorithm_parameters_yaml=algorithm_parameters_yaml,
            platform_id=platform_id,
            source_data_software_disp_version=source_data_software_disp_version,
            cal_disp_software_version=cal_disp_software_version,
            venti_software_version=venti_software_version,
            product_pixel_coordinate_convention=product_pixel_coordinate_convention,
            ceos_atmospheric_phase_correction=ceos_atmospheric_phase_correction,
            ceos_gridding_convention=ceos_gridding_convention,
            ceos_product_measurement_projection=ceos_product_measurement_projection,
            ceos_ionospheric_phase_correction=ceos_ionospheric_phase_correction,
            ceos_noise_removal=ceos_noise_removal,
            algorithm_theoretical_basis_document_doi=algorithm_theoretical_basis_document_doi,
            product_landing_page_doi=product_landing_page_doi,
            product_specification_document_id=product_specification_document_id,
            pge_runconfig=pge_runconfig,
        )
        ds_meta.to_netcdf(self.path, mode="a", group="metadata", engine="h5netcdf")

    def add_auxiliary(
        self,
        model_3d: dict[str, xr.DataArray] | None = None,
        model_3d_std: dict[str, xr.DataArray] | None = None,
        spatial_ref: xr.DataArray | None = None,
    ) -> None:
        """Add auxiliary group with 3D displacement model to existing product.

        Parameters
        ----------
        model_3d : dict[str, xr.DataArray] or None, optional
            3D displacement components with keys:
            "north_south", "east_west", "up_down".
        model_3d_std : dict[str, xr.DataArray] or None, optional
            3D displacement uncertainties.
        spatial_ref : xr.DataArray or None, optional
            Spatial reference.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.
        ValueError
            If neither model_3d nor model_3d_std provided.

        Examples
        --------
        >>> cal.add_auxiliary(
        ...     model_3d={
        ...         "north_south": ns_array,
        ...         "east_west": ew_array,
        ...         "up_down": ud_array,
        ...     }
        ... )

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        if not model_3d and not model_3d_std:
            raise ValueError("Must provide at least one of model_3d or model_3d_std")

        ds_aux = build_auxiliary_dataset(
            model_3d=model_3d,
            model_3d_std=model_3d_std,
            spatial_ref=spatial_ref,
        )
        ds_aux.to_netcdf(self.path, mode="a", group="auxiliary", engine="h5netcdf")

    def open_dataset(self, group: str | None = None) -> xr.Dataset:
        """Open calibration dataset.

        Parameters
        ----------
        group : str or None, optional
            Group to open: None for main, "auxiliary" for 3D model.
            Default is None (main group).

        Returns
        -------
        xr.Dataset
            Dataset containing requested group.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.

        Examples
        --------
        >>> # Open main calibration (full resolution)
        >>> ds_main = cal.open_dataset()
        >>> calibration = ds_main["calibration"]

        >>> # Open auxiliary group (coarse resolution)
        >>> ds_model = cal.open_dataset(group="auxiliary")
        >>> model_up = ds_model["up_down"]

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        if group == "auxiliary":
            return xr.open_dataset(self.path, group="auxiliary", engine="h5netcdf")

        return xr.open_dataset(self.path, engine="h5netcdf")

    def open_auxiliary(self) -> xr.Dataset:
        """Open auxiliary group dataset.

        Returns
        -------
        xr.Dataset
            Dataset containing 3D displacement model at coarse resolution.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.
        ValueError
            If auxiliary group does not exist.

        Examples
        --------
        >>> ds_model = cal.open_auxiliary()
        >>> disp_ns = ds_model["north_south"]
        >>> disp_ew = ds_model["east_west"]
        >>> disp_up = ds_model["up_down"]

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        try:
            return xr.open_dataset(self.path, group="auxiliary", engine="h5netcdf")
        except (OSError, ValueError) as e:
            raise ValueError(
                f"auxiliary group not found in {self.filename}. "
                "Product may not contain 3D displacement model."
            ) from e

    def open_identification(self) -> xr.Dataset:
        """Open identification group dataset.

        Returns
        -------
        xr.Dataset
            Dataset containing product identification metadata.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.
        ValueError
            If identification group does not exist.

        Examples
        --------
        >>> ds_id = cal.open_identification()
        >>> frame_id = ds_id["frame_id"].item()
        >>> pge_config = ds_id["pge_runconfig"].item()

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        try:
            return xr.open_dataset(self.path, group="identification", engine="h5netcdf")
        except (OSError, ValueError) as e:
            raise ValueError(
                f"identification group not found in {self.filename}."
            ) from e

    def open_metadata(self) -> xr.Dataset:
        """Open metadata group dataset.

        Returns
        -------
        xr.Dataset
            Dataset containing processing and algorithm metadata.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.
        ValueError
            If metadata group does not exist.

        Examples
        --------
        >>> ds_meta = cal.open_metadata()
        >>> platform = ds_meta["platform_id"].item()

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        try:
            return xr.open_dataset(self.path, group="metadata", engine="h5netcdf")
        except (OSError, ValueError) as e:
            raise ValueError(f"metadata group not found in {self.filename}.") from e

    def has_auxiliary(self) -> bool:
        """Check if product contains auxiliary group.

        Returns
        -------
        bool
            True if auxiliary group exists.

        """
        try:
            self.open_auxiliary()
            return True
        except (FileNotFoundError, ValueError):
            return False

    def has_metadata(self) -> bool:
        """Check if product contains metadata group."""
        try:
            self.open_metadata()
            return True
        except (FileNotFoundError, ValueError):
            return False

    def has_identification(self) -> bool:
        """Check if product contains identification group."""
        try:
            self.open_identification()
            return True
        except (FileNotFoundError, ValueError):
            return False

    def get_epsg(self) -> int | None:
        """Get EPSG code from spatial reference.

        Returns
        -------
        int or None
            EPSG code if available, None otherwise.

        """
        ds = self.open_dataset()

        if "spatial_ref" in ds:
            crs_wkt = ds.spatial_ref.attrs.get("crs_wkt")
            if crs_wkt:
                crs = CRS.from_wkt(crs_wkt)
                return crs.to_epsg()

        return None

    def get_bounds(self) -> dict[str, float]:
        """Get bounds in native projection.

        Returns
        -------
        dict[str, float]
            Bounds with keys: left, bottom, right, top.

        """
        ds = self.open_dataset()

        x = ds.x.values
        y = ds.y.values

        return {
            "left": float(x.min()),
            "bottom": float(y.min()),
            "right": float(x.max()),
            "top": float(y.max()),
        }

    def get_bounds_wgs84(self) -> dict[str, float]:
        """Get bounds transformed to WGS84.

        Returns
        -------
        dict[str, float]
            Bounds in WGS84 with keys: west, south, east, north.

        Raises
        ------
        ValueError
            If spatial_ref or crs_wkt is missing.

        """
        ds = self.open_dataset()

        x = ds.x.values
        y = ds.y.values
        left = float(x.min())
        bottom = float(y.min())
        right = float(x.max())
        top = float(y.max())

        if "spatial_ref" not in ds:
            raise ValueError("Dataset missing spatial_ref")

        crs_wkt = ds.spatial_ref.attrs.get("crs_wkt")
        if not crs_wkt:
            raise ValueError("spatial_ref missing crs_wkt")

        src_crs = CRS.from_wkt(crs_wkt)

        west, south, east, north = transform_bounds(
            src_crs,
            CRS.from_epsg(4326),
            left,
            bottom,
            right,
            top,
        )

        return {
            "west": west,
            "south": south,
            "east": east,
            "north": north,
        }

    def to_geotiff(
        self,
        layer: str,
        output_path: Path | str,
        time_index: int = 0,
        group: str | None = None,
        compress: str = "DEFLATE",
        **kwargs,
    ) -> Path:
        """Export layer to GeoTIFF.

        Parameters
        ----------
        layer : str
            Name of layer to export.
        output_path : Path or str
            Output GeoTIFF path.
        time_index : int, optional
            Time index to export (for 3D data). Default is 0 (first time).
        group : str or None, optional
            Group containing layer. Default is None (main group).
        compress : str, optional
            Compression method. Default is "DEFLATE".
        **kwargs
            Additional rasterio creation options.

        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        ds = self.open_dataset(group=group)

        if layer not in ds:
            available = list(ds.data_vars)
            group_str = f" in {group} group" if group else ""
            raise ValueError(
                f"Layer '{layer}' not found{group_str}. Available: {available}"
            )

        da = ds[layer]

        # Handle 3D data (time, y, x) by selecting time slice
        if "time" in da.dims:
            if time_index >= len(da.time):
                raise ValueError(
                    f"time_index {time_index} out of range for {len(da.time)} time"
                    " steps"
                )
            da = da.isel(time=time_index)

        data = da.values

        # Extract spatial information using imported helpers
        if "spatial_ref" in ds:
            transform = get_transform(ds)
            crs = get_crs(ds)
        else:
            transform = compute_transform_from_coords(ds.x.values, ds.y.values)
            crs = None

        # Write GeoTIFF
        profile = {
            "driver": "GTiff",
            "height": data.shape[0],
            "width": data.shape[1],
            "count": 1,
            "dtype": np.float32,
            "transform": transform,
            "compress": compress,
            "tiled": True,
            "blockxsize": 512,
            "blockysize": 512,
            **kwargs,
        }

        if crs:
            profile["crs"] = crs

        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(data.astype(np.float32), 1)
            dst.set_band_description(1, layer)

            # Add OPERA metadata tags
            dst.update_tags(
                product_type=f"OPERA_L4_CAL-DISP-{self.sensor}",
                sensor=self.sensor,
                frame_id=self.frame_id,
                polarization=self.polarization,
                reference_date=self.reference_date.isoformat(),
                secondary_date=self.secondary_date.isoformat(),
                layer=layer,
                group=group if group else "main",
            )

        return output_path

    def open_all_groups(self) -> dict[str, xr.Dataset]:
        """Open all available groups.

        Returns
        -------
        dict[str, xr.Dataset]
            Dictionary with keys "main", and optionally "identification",
            "metadata", and "auxiliary".

        """
        groups = {"main": self.open_dataset()}

        if self.has_identification():
            groups["identification"] = self.open_identification()

        if self.has_metadata():
            groups["metadata"] = self.open_metadata()

        if self.has_auxiliary():
            groups["auxiliary"] = self.open_auxiliary()

        return groups

    @property
    def baseline_days(self) -> int:
        """Temporal baseline in days."""
        return (self.secondary_date - self.reference_date).days

    @property
    def filename(self) -> str:
        """Product filename."""
        return self.path.name

    @property
    def exists(self) -> bool:
        """Check if product file exists."""
        return self.path.exists()

    def __repr__(self) -> str:
        """Return a string representation."""
        model_str = "+auxiliary" if self.exists and self.has_auxiliary() else ""
        return (
            f"CalProduct(sensor={self.sensor}, frame={self.frame_id}, "
            f"{self.reference_date.date()} → {self.secondary_date.date()}, "
            f"{self.polarization}{model_str})"
        )
