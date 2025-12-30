import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Literal

import rasterio
import xarray as xr
from rasterio.crs import CRS
from rasterio.transform import Affine
from rasterio.warp import transform_bounds


@dataclass
class DispProduct:
    """OPERA DISP-S1 displacement product.

    Represents a Level-3 interferometric displacement product from
    the OPERA DISP-S1 archive. Products contain displacement measurements
    between two Sentinel-1 acquisition dates.

    Parameters
    ----------
    path : Path
        Path to the NetCDF product file.
    frame_id : int
        OPERA frame identifier (e.g., 8882).
    primary_date : datetime
        Earlier acquisition date (reference).
    secondary_date : datetime
        Later acquisition date.
    polarization : str
        Radar polarization (e.g., "VV", "VH").
    version : str
        Product version string (e.g., "1.0").
    production_date : datetime
        Date when product was generated.
    mode : str, optional
        Acquisition mode (e.g., "IW"). Default is "IW".

    Examples
    --------
    >>> path = Path(
    "OPERA_L3_DISP-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0_20251027T005420Z.nc")
    >>> product = DispProduct.from_path(path)
    >>> product.frame_id
    8882

    # Get reference point
    >>> row, col = product.get_reference_point_index()
    >>> lat, lon = product.get_reference_point_latlon()
    >>> print(f"Reference at ({row}, {col}): {lat:.4f}°N, {lon:.4f}°E")

    """

    path: Path
    frame_id: int
    primary_date: datetime
    secondary_date: datetime
    polarization: str
    version: str
    production_date: datetime
    mode: str = "IW"

    # Filename pattern for OPERA DISP-S1 products
    _PATTERN = re.compile(
        r"OPERA_L3_DISP-S1_"
        r"(?P<mode>\w+)_"
        r"F(?P<frame_id>\d+)_"
        r"(?P<pol>\w+)_"
        r"(?P<primary>\d{8}T\d{6}Z)_"
        r"(?P<secondary>\d{8}T\d{6}Z)_"
        r"v(?P<version>[\d.]+)_"
        r"(?P<production>\d{8}T\d{6}Z)"
        r"\.nc$"
    )

    # Main dataset layers
    DISPLACEMENT_LAYERS = [
        "displacement",
        "short_wavelength_displacement",
        "recommended_mask",
        "connected_component_labels",
        "temporal_coherence",
        "estimated_phase_quality",
        "persistent_scatterer_mask",
        "shp_counts",
        "water_mask",
        "phase_similarity",
        "timeseries_inversion_residuals",
    ]

    # Corrections group layers
    CORRECTION_LAYERS = [
        "ionospheric_delay",
        "solid_earth_tide",
        "perpendicular_baseline",
    ]

    def __post_init__(self) -> None:
        """Validate product after construction."""
        self.path = Path(self.path)

        if self.frame_id <= 0:
            raise ValueError(f"frame_id must be positive, got {self.frame_id}")

        if self.secondary_date <= self.primary_date:
            raise ValueError(
                f"Secondary date ({self.secondary_date}) must be after "
                f"primary date ({self.primary_date})"
            )

        if self.polarization not in {"VV", "VH", "HH", "HV"}:
            raise ValueError(f"Invalid polarization: {self.polarization}")

    @classmethod
    def from_path(cls, path: Path | str) -> "DispProduct":
        """Parse product metadata from filename.

        Parameters
        ----------
        path : Path or str
            Path to OPERA DISP-S1 NetCDF file.

        Returns
        -------
        DispProduct
            Parsed product instance.

        Raises
        ------
        ValueError
            If filename doesn't match expected OPERA DISP-S1 format.

        """
        path = Path(path)
        match = cls._PATTERN.match(path.name)

        if not match:
            raise ValueError(
                f"Filename does not match OPERA DISP-S1 pattern: {path.name}"
            )

        return cls(
            path=path,
            frame_id=int(match.group("frame_id")),
            primary_date=datetime.strptime(match.group("primary"), "%Y%m%dT%H%M%SZ"),
            secondary_date=datetime.strptime(
                match.group("secondary"), "%Y%m%dT%H%M%SZ"
            ),
            polarization=match.group("pol"),
            version=match.group("version"),
            production_date=datetime.strptime(
                match.group("production"), "%Y%m%dT%H%M%SZ"
            ),
            mode=match.group("mode"),
        )

    def open_dataset(
        self, group: Literal["main", "corrections"] | None = None
    ) -> xr.Dataset:
        """Open dataset.

        Parameters
        ----------
        group : {"main", "corrections"} or None, optional
            Which group to open. If None, opens main group. Default is None.

        Returns
        -------
        xr.Dataset
            Dataset containing displacement and quality layers.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        if group == "corrections":
            return xr.open_dataset(self.path, group="corrections", engine="h5netcdf")

        return xr.open_dataset(self.path, engine="h5netcdf")

    def open_corrections(self) -> xr.Dataset:
        """Open corrections group dataset.

        Returns
        -------
        xr.Dataset
            Corrections dataset containing ionospheric delay, solid earth tide, etc.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.

        """
        return self.open_dataset(group="corrections")

    def get_epsg(self) -> int | None:
        """Get EPSG code from spatial reference.

        Returns
        -------
        int or None
            EPSG code if found, None otherwise.

        Examples
        --------
        >>> product.get_epsg()
        32615

        """
        ds = self.open_dataset()
        crs = self._get_crs(ds)

        # Parse EPSG from CRS
        raster_crs = CRS.from_wkt(crs)
        return raster_crs.to_epsg()

    def get_bounds(self) -> dict[str, float]:
        """Get bounds in native projection coordinates.

        Returns
        -------
        dict[str, float]
            Dictionary with keys: left, bottom, right, top.

        Examples
        --------
        >>> bounds = product.get_bounds()
        >>> bounds
        {'left': 71970.0, 'bottom': 3153930.0, 'right': 355890.0, 'top': 3385920.0}

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
        """Get bounds transformed to WGS84 (EPSG:4326).

        Returns
        -------
        dict[str, float]
            Dictionary with keys: west, south, east, north in decimal degrees.

        Examples
        --------
        >>> bounds = product.get_bounds_wgs84()
        >>> bounds
        {'west': -95.567, 'south': 28.486, 'east': -93.212, 'north': 30.845}

        """
        ds = self.open_dataset()

        # Get native bounds
        x = ds.x.values
        y = ds.y.values
        left = float(x.min())
        bottom = float(y.min())
        right = float(x.max())
        top = float(y.max())

        # Get native CRS
        crs_wkt = self._get_crs(ds)
        src_crs = CRS.from_wkt(crs_wkt)

        # Transform to WGS84
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

    def get_reference_point_index(self) -> tuple[int, int]:
        """Get reference point pixel indices.

        The reference point is where the phase was set to zero during
        processing. This is stored in the corrections group.

        Returns
        -------
        tuple[int, int]
            Row and column indices (row, col) of reference point.

        Raises
        ------
        ValueError
            If reference_point variable not found or missing attributes.

        Examples
        --------
        >>> row, col = product.get_reference_point_index()
        >>> print(f"Reference point at pixel ({row}, {col})")

        """
        ds = self.open_corrections()

        if "reference_point" not in ds:
            raise ValueError("reference_point variable not found in corrections group")

        ref_attrs = ds.reference_point.attrs

        if "rows" not in ref_attrs or "cols" not in ref_attrs:
            raise ValueError("reference_point missing 'rows' or 'cols' attributes")

        row = int(ref_attrs["rows"])
        col = int(ref_attrs["cols"])

        return (row, col)

    def get_reference_point_latlon(self) -> tuple[float, float]:
        """Get reference point geographic coordinates.

        Returns latitude and longitude of the reference point in WGS84.

        Returns
        -------
        tuple[float, float]
            Latitude and longitude in decimal degrees (lat, lon).

        Raises
        ------
        ValueError
            If reference_point variable not found or missing attributes.

        Examples
        --------
        >>> lat, lon = product.get_reference_point_latlon()
        >>> print(f"Reference point: {lat:.6f}°N, {lon:.6f}°E")

        """
        ds = self.open_corrections()

        if "reference_point" not in ds:
            raise ValueError("reference_point variable not found in corrections group")

        ref_attrs = ds.reference_point.attrs

        if "latitudes" not in ref_attrs or "longitudes" not in ref_attrs:
            raise ValueError(
                "reference_point missing 'latitudes' or 'longitudes' attributes"
            )

        lat = float(ref_attrs["latitudes"])
        lon = float(ref_attrs["longitudes"])

        return (lat, lon)

    def to_geotiff(
        self,
        layer: str,
        output_path: Path | str,
        group: Literal["main", "corrections"] = "main",
        compress: str = "DEFLATE",
        **kwargs,
    ) -> Path:
        """Export layer to optimized GeoTIFF.

        Parameters
        ----------
        layer : str
            Name of layer to export (e.g., "displacement", "ionospheric_delay").
        output_path : Path or str
            Output GeoTIFF path.
        group : {"main", "corrections"}, optional
            Which group to read from. Default is "main".
        compress : str, optional
            Compression method. Default is "DEFLATE".
        **kwargs
            Additional rasterio creation options.

        Returns
        -------
        Path
            Path to created GeoTIFF.

        Raises
        ------
        ValueError
            If layer not found in specified group.

        Examples
        --------
        >>> product.to_geotiff("displacement", "disp.tif")
        >>> product.to_geotiff("ionospheric_delay", "iono.tif", group="corrections")

        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Open appropriate dataset
        ds = self.open_dataset(group=group if group == "corrections" else None)

        if layer not in ds:
            available = list(ds.data_vars)
            raise ValueError(
                f"Layer '{layer}' not found in {group} group. "
                f"Available layers: {available}"
            )

        # Get data array
        da = ds[layer]

        # Extract transform from spatial_ref
        transform = self._get_transform(ds)

        # Get CRS
        crs = self._get_crs(ds)

        # Prepare data - handle (y, x) or (time, y, x) shapes
        if da.ndim == 3:
            # Take first time slice if 3D
            data = da.isel(time=0).values
        else:
            data = da.values

        # Write GeoTIFF
        profile = {
            "driver": "GTiff",
            "height": data.shape[0],
            "width": data.shape[1],
            "count": 1,
            "dtype": data.dtype,
            "crs": crs,
            "transform": transform,
            "compress": compress,
            "tiled": True,
            "blockxsize": 512,
            "blockysize": 512,
            **kwargs,
        }

        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(data, 1)
            dst.set_band_description(1, layer)

        return output_path

    def _get_transform(self, ds: xr.Dataset) -> Affine:
        """Extract affine transform from dataset."""
        gt = ds.spatial_ref.attrs.get("GeoTransform")
        if gt is None:
            raise ValueError("No GeoTransform found in spatial_ref")

        # Parse string like "71970.0 30.0 0.0 3385920.0 0.0 -30.0"
        vals = [float(x) for x in gt.split()]
        return Affine(vals[1], vals[2], vals[0], vals[4], vals[5], vals[3])

    def _get_crs(self, ds: xr.Dataset) -> str:
        """Extract CRS from dataset."""
        crs_wkt = ds.spatial_ref.attrs.get("crs_wkt")
        if crs_wkt is None:
            raise ValueError("No crs_wkt found in spatial_ref")
        return crs_wkt

    @property
    def baseline_days(self) -> int:
        """Temporal baseline in days between acquisitions."""
        return (self.secondary_date - self.primary_date).days

    @property
    def filename(self) -> str:
        """Product filename."""
        return self.path.name

    @property
    def exists(self) -> bool:
        """Check if product file exists on disk."""
        return self.path.exists()

    def __repr__(self) -> str:
        """Concise string representation."""
        return (
            f"DispProduct(frame={self.frame_id}, "
            f"{self.primary_date.date()} → {self.secondary_date.date()}, "
            f"{self.polarization})"
        )
