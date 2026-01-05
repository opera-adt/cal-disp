from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
import rasterio
import xarray as xr
from rasterio.crs import CRS
from rasterio.transform import Affine
from rasterio.warp import transform_bounds


@dataclass
class StaticLayer:
    """OPERA DISP-S1-STATIC layer.

    Represents a single static layer (DEM, incidence angle, LOS vectors, etc.)
    used as input for DISP-S1 processing. These are frame-specific GeoTIFF
    files that don't change over time.

    Examples
    --------
    >>> path = Path("OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_dem.tif")
    >>> layer = StaticLayer.from_path(path)
    >>> layer.frame_id
    8882

    # Read LOS components
    >>> los_layer = StaticLayer.from_path("..._line_of_sight_enu.tif")
    >>> bands = los_layer.read_bands()
    >>> east, north, up = bands[0], bands[1], bands[2]

    """

    path: Path
    frame_id: int
    reference_date: datetime
    satellite: str
    version: str
    layer_type: str

    _PATTERN = re.compile(
        r"OPERA_L3_DISP-S1-STATIC_"
        r"F(?P<frame_id>\d+)_"
        r"(?P<date>\d{8})_"
        r"(?P<satellite>S1[AB])_"
        r"v(?P<version>[\d.]+)_"
        r"(?P<layer>[\w_]+)"
        r"\.tif$"
    )

    LAYER_TYPES = [
        "dem",
        "line_of_sight_enu",
        "layover_shadow_mask",
    ]

    def __post_init__(self) -> None:
        """Validate layer after construction."""
        self.path = Path(self.path)

        if self.frame_id <= 0:
            raise ValueError(f"frame_id must be positive, got {self.frame_id}")

    @classmethod
    def from_path(cls, path: Path | str) -> "StaticLayer":
        """Parse layer metadata from filename."""
        path = Path(path)
        match = cls._PATTERN.match(path.name)

        if not match:
            raise ValueError(
                f"Filename does not match OPERA DISP-S1-STATIC pattern: {path.name}"
            )

        return cls(
            path=path,
            frame_id=int(match.group("frame_id")),
            reference_date=datetime.strptime(match.group("date"), "%Y%m%d"),
            satellite=match.group("satellite"),
            version=match.group("version"),
            layer_type=match.group("layer"),
        )

    @property
    def num_bands(self) -> int:
        """Get number of bands in layer."""
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            return src.count

    def read(self, band: int = 1, masked: bool = True) -> np.ndarray:
        """Read single band data."""
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            if band < 1 or band > src.count:
                raise ValueError(
                    f"Band {band} out of range. File has {src.count} bands."
                )

            data = src.read(band)

            if masked and src.nodata is not None:
                data = np.ma.masked_equal(data, src.nodata)

        return data

    def read_bands(self, masked: bool = True) -> list[np.ndarray]:
        """Read all bands.

        Parameters
        ----------
        masked : bool, optional
            If True, return masked arrays with nodata values masked.
            Default is True.

        Returns
        -------
        list[np.ndarray]
            List of arrays, one per band.

        Examples
        --------
        >>> # Read DEM (single band)
        >>> dem_layer = StaticLayer.from_path("..._dem.tif")
        >>> bands = dem_layer.read_bands()
        >>> dem = bands[0]

        >>> # Read LOS vectors (three bands)
        >>> los_layer = StaticLayer.from_path("..._line_of_sight_enu.tif")
        >>> bands = los_layer.read_bands()
        >>> east, north, up = bands[0], bands[1], bands[2]

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            bands = []
            for band_idx in range(1, src.count + 1):
                data = src.read(band_idx)
                if masked and src.nodata is not None:
                    data = np.ma.masked_equal(data, src.nodata)
                bands.append(data)

        return bands

    def to_dataset(self) -> xr.Dataset:
        """Convert raster to xarray Dataset."""
        # Get shape and transform using existing methods
        height, width = self.get_shape()
        transform = self.get_transform()

        # Generate x and y coordinates from transform
        x_coords = np.arange(width) * transform[0] + transform[2] + transform[0] / 2
        y_coords = np.arange(height) * transform[4] + transform[5] + transform[4] / 2

        # Create coordinates dict
        coords = {
            "y": (["y"], y_coords),
            "x": (["x"], x_coords),
        }

        # Create data variables
        data_vars = {}
        nodata = self.get_nodata()

        if self.num_bands == 1:
            # Single band - use layer_type as variable name
            data = self.read(band=1, masked=False)
            if nodata is not None:
                data = np.where(data == nodata, np.nan, data)
            data_vars[self.layer_type] = (["y", "x"], data)

        elif self.layer_type == "line_of_sight_enu" and self.num_bands == 3:
            # Special case: LOS vectors - create three variables
            band_names = ["los_east", "los_north", "los_up"]
            bands = self.read_bands(masked=False)

            for name, data in zip(band_names, bands):
                if nodata is not None:
                    data = np.where(data == nodata, np.nan, data)
                data_vars[name] = (["y", "x"], data)

        else:
            # Generic multi-band: use band dimension
            bands = self.read_bands(masked=False)
            all_data = np.stack(bands)

            if nodata is not None:
                all_data = np.where(all_data == nodata, np.nan, all_data)

            coords["band"] = (["band"], np.arange(1, self.num_bands + 1))
            data_vars[self.layer_type] = (["band", "y", "x"], all_data)

        # Create dataset
        ds = xr.Dataset(data_vars=data_vars, coords=coords)

        # Add attributes
        ds.attrs["frame_id"] = self.frame_id
        ds.attrs["satellite"] = self.satellite
        ds.attrs["version"] = self.version
        ds.attrs["reference_date"] = self.reference_date.isoformat()
        ds.attrs["layer_type"] = self.layer_type

        # Add CRS information using existing method
        crs = self.get_crs()
        if crs:
            ds.attrs["crs_wkt"] = crs.to_wkt()
            epsg = self.get_epsg()
            if epsg:
                ds.attrs["epsg"] = epsg

        # Add transform
        ds.attrs["transform"] = list(transform)
        ds.attrs["nodata"] = nodata

        # Add coordinate attributes
        ds["x"].attrs["standard_name"] = "projection_x_coordinate"
        ds["x"].attrs["long_name"] = "x coordinate of projection"
        ds["x"].attrs["units"] = "m"

        ds["y"].attrs["standard_name"] = "projection_y_coordinate"
        ds["y"].attrs["long_name"] = "y coordinate of projection"
        ds["y"].attrs["units"] = "m"

        # Add variable-specific attributes
        if self.layer_type == "dem":
            ds["dem"].attrs["units"] = "m"
            ds["dem"].attrs["long_name"] = "Digital Elevation Model"
        elif self.layer_type == "line_of_sight_enu":
            ds["los_east"].attrs["long_name"] = "LOS unit vector - East component"
            ds["los_north"].attrs["long_name"] = "LOS unit vector - North component"
            ds["los_up"].attrs["long_name"] = "LOS unit vector - Up component"
        elif self.layer_type == "layover_shadow_mask":
            ds["layover_shadow_mask"].attrs["long_name"] = "Layover and Shadow Mask"
            ds["layover_shadow_mask"].attrs[
                "description"
            ] = "0=valid, 1=layover, 2=shadow"

        return ds

    def compute_incidence_angle(
        self,
        fill_value: float = 0.0,
        dtype: np.dtype = np.float32,
    ) -> np.ndarray:
        """Compute incidence angle from LOS up component.

        Only valid for line_of_sight_enu layers.

        Parameters
        ----------
        fill_value : float, optional
            Value to use for masked/invalid pixels. Default is 0.0.
        dtype : np.dtype, optional
            Output data type. Default is np.float32.

        Returns
        -------
        np.ndarray
            Incidence angle in degrees (0-90°).

        Raises
        ------
        ValueError
            If not a line_of_sight_enu layer or wrong number of bands.

        """
        if self.layer_type != "line_of_sight_enu":
            raise ValueError(
                "compute_incidence_angle() only valid for line_of_sight_enu layers, "
                f"got {self.layer_type}"
            )

        if self.num_bands != 3:
            raise ValueError(f"Expected 3 bands for LOS ENU, got {self.num_bands}")

        # Get LOS up component (band 3)
        bands = self.read_bands()
        los_up = bands[2]

        # Handle masked arrays
        if isinstance(los_up, np.ma.MaskedArray):
            los_up_data = los_up.data
            mask = los_up.mask
        else:
            los_up_data = los_up
            mask = None

        # Clip to valid range [-1, 1] to avoid arccos domain errors
        los_up_clipped = np.clip(los_up_data, -1.0, 1.0)

        # Compute incidence angle in degrees
        incidence_angle = np.rad2deg(np.arccos(los_up_clipped))

        # Apply mask if present
        if mask is not None:
            incidence_angle = np.where(mask, fill_value, incidence_angle)

        return incidence_angle.astype(dtype)

    def export_incidence_angle(
        self,
        output_path: Path | str,
        fill_value: float = 0.0,
        nodata: float | None = 0.0,
        compress: str = "DEFLATE",
        **kwargs,
    ) -> Path:
        """Compute and export incidence angle to GeoTIFF."""
        if self.layer_type != "line_of_sight_enu":
            raise ValueError(
                "export_incidence_angle() only valid for line_of_sight_enu layers, "
                f"got {self.layer_type}"
            )

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Compute incidence angle
        incidence_angle = self.compute_incidence_angle(fill_value=fill_value)

        # Write GeoTIFF
        profile = {
            "driver": "GTiff",
            "height": incidence_angle.shape[0],
            "width": incidence_angle.shape[1],
            "count": 1,
            "dtype": incidence_angle.dtype,
            "crs": self.get_crs(),
            "transform": self.get_transform(),
            "nodata": nodata,
            "compress": compress,
            "tiled": True,
            "blockxsize": 512,
            "blockysize": 512,
            **kwargs,
        }

        with rasterio.open(output_path, "w", **profile) as dst:
            dst.write(incidence_angle, 1)
            dst.set_band_description(1, "Incidence angle (degrees)")
            dst.update_tags(
                1,
                units="degrees",
                description="Incidence angle computed from LOS up component",
            )

        return output_path

    def get_profile(self) -> dict:
        """Get rasterio profile."""
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            return src.profile

    def get_transform(self) -> Affine:
        """Get affine transform."""
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            return src.transform

    def get_crs(self) -> CRS:
        """Get coordinate reference system."""
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            return src.crs

    def get_epsg(self) -> int | None:
        """Get EPSG code."""
        crs = self.get_crs()
        return crs.to_epsg() if crs else None

    def get_bounds(self) -> dict[str, float]:
        """Get bounds in native projection."""
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            bounds = src.bounds
            return {
                "left": bounds.left,
                "bottom": bounds.bottom,
                "right": bounds.right,
                "top": bounds.top,
            }

    def get_bounds_wgs84(self) -> dict[str, float]:
        """Get bounds transformed to WGS84."""
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            bounds = src.bounds
            west, south, east, north = transform_bounds(
                src.crs,
                CRS.from_epsg(4326),
                bounds.left,
                bounds.bottom,
                bounds.right,
                bounds.top,
            )

        return {
            "west": west,
            "south": south,
            "east": east,
            "north": north,
        }

    def get_shape(self) -> tuple[int, int]:
        """Get array shape."""
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            return (src.height, src.width)

    def get_nodata(self) -> float | None:
        """Get nodata value."""
        if not self.path.exists():
            raise FileNotFoundError(f"Layer file not found: {self.path}")

        with rasterio.open(self.path) as src:
            return src.nodata

    @property
    def filename(self) -> str:
        """Layer filename."""
        return self.path.name

    @property
    def exists(self) -> bool:
        """Check if layer file exists."""
        return self.path.exists()

    def __repr__(self) -> str:
        """Return a string representation."""
        band_info = f", bands={self.num_bands}" if self.exists else ""
        return f"StaticLayer(frame={self.frame_id}, layer={self.layer_type}{band_info})"
