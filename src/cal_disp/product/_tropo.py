from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import xarray as xr
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.warp import transform_bounds as rio_transform_bounds
from scipy.interpolate import RegularGridInterpolator


@dataclass
class TropoProduct:
    """OPERA TROPO-ZENITH tropospheric delay product.

    Minimal class for managing OPERA tropospheric products.
    Processing functions are standalone for composability.

    Examples
    --------
    >>> product = TropoProduct.from_path("tropo.nc")
    >>> ds = product.open_dataset()
    >>> total = product.get_total_delay()

    """

    path: Path
    date: datetime
    production_date: datetime
    model: str
    version: str

    _PATTERN = re.compile(
        r"OPERA_L4_TROPO-ZENITH_"
        r"(?P<date>\d{8}T\d{6}Z)_"
        r"(?P<production>\d{8}T\d{6}Z)_"
        r"(?P<model>\w+)_"
        r"v(?P<version>[\d.]+)"
        r"\.nc$"
    )

    TROPO_LAYERS = [
        "wet_delay",
        "hydrostatic_delay",
    ]

    def __post_init__(self) -> None:
        """Validate product after construction."""
        self.path = Path(self.path)

        if self.production_date < self.date:
            raise ValueError(
                f"Production date ({self.production_date}) cannot be before "
                f"model date ({self.date})"
            )

    @classmethod
    def from_path(cls, path: Path | str) -> "TropoProduct":
        """Parse product metadata from filename."""
        path = Path(path)
        match = cls._PATTERN.match(path.name)

        if not match:
            raise ValueError(
                f"Filename does not match OPERA TROPO-ZENITH pattern: {path.name}"
            )

        return cls(
            path=path,
            date=datetime.strptime(match.group("date"), "%Y%m%dT%H%M%SZ"),
            production_date=datetime.strptime(
                match.group("production"), "%Y%m%dT%H%M%SZ"
            ),
            model=match.group("model"),
            version=match.group("version"),
        )

    def matches_date(self, target_date: datetime, hours: float = 6.0) -> bool:
        """Check if product date is within time window of target date."""
        delta = abs(self.date - target_date)
        return delta <= timedelta(hours=hours)

    def open_dataset(
        self,
        bounds: tuple[float, float, float, float] | None = None,
        max_height: float | None = None,
        bounds_crs: str = "EPSG:4326",
        bounds_buffer: float = 0.0,
    ) -> xr.Dataset:
        """Open tropospheric delay dataset with optional subsetting.

        Parameters
        ----------
        bounds : tuple[float, float, float, float] or None, optional
            Spatial bounds as (west, south, east, north). Default is None.
        max_height : float or None, optional
            Maximum height in meters. Default is None.
        bounds_crs : str, optional
            CRS of bounds. Default is "EPSG:4326".
        bounds_buffer : float, optional
            Buffer to add to bounds in degrees (for lat/lon) or meters
            (for projected CRS). Default is 0.0. Useful value: 0.2 for lat/lon.

        Returns
        -------
        xr.Dataset
            Tropospheric delay dataset.

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        ds = xr.open_dataset(self.path, engine="h5netcdf")

        if bounds is not None:
            # Apply buffer if requested
            if bounds_buffer > 0:
                west, south, east, north = bounds
                bounds = (
                    west - bounds_buffer,
                    south - bounds_buffer,
                    east + bounds_buffer,
                    north + bounds_buffer,
                )
            ds = self._subset_spatial(ds, bounds, bounds_crs)

        if max_height is not None:
            ds = self._subset_height(ds, max_height)

        return ds

    def get_total_delay(
        self,
        time_idx: int = 0,
        bounds: tuple[float, float, float, float] | None = None,
        max_height: float | None = None,
        bounds_crs: str = "EPSG:4326",
        bounds_buffer: float = 0.0,
    ) -> xr.DataArray:
        """Get total zenith delay (wet + hydrostatic).

        Computes total delay as sum of wet and hydrostatic components.

        Parameters
        ----------
        time_idx : int, optional
            Time index to extract. Default is 0.
        bounds : tuple[float, float, float, float] or None, optional
            Spatial bounds for subsetting. Default is None.
        max_height : float or None, optional
            Maximum height for subsetting. Default is None.
        bounds_crs : str, optional
            CRS of bounds. Default is "EPSG:4326".
        bounds_buffer : float, optional
            Buffer to add to bounds. Default is 0.0.

        Returns
        -------
        xr.DataArray
            Total zenith delay with dimensions (height, latitude, longitude).

        """
        ds = self.open_dataset(
            bounds=bounds,
            max_height=max_height,
            bounds_crs=bounds_crs,
            bounds_buffer=bounds_buffer,
        )

        # Compute total delay from wet + hydrostatic
        if "wet_delay" not in ds:
            raise ValueError("wet_delay not found in dataset")
        if "hydrostatic_delay" not in ds:
            raise ValueError("hydrostatic_delay not found in dataset")

        da = ds["wet_delay"] + ds["hydrostatic_delay"]
        da.name = "zenith_total_delay"
        da.attrs.update(
            {
                "long_name": "Total zenith tropospheric delay",
                "units": "meters",
                "description": "Sum of wet and hydrostatic components",
            }
        )

        # Preserve spatial reference
        if "spatial_ref" in ds:
            da = da.assign_coords({"spatial_ref": ds["spatial_ref"]})

        # Handle time dimension
        if "time" in da.dims:
            n_times = len(da.time)
            if abs(time_idx) >= n_times:
                raise ValueError(
                    f"time_idx {time_idx} out of range for {n_times} timesteps"
                )
            da = da.isel(time=time_idx)

        return da

    def _subset_spatial(
        self,
        ds: xr.Dataset,
        bounds: tuple[float, float, float, float],
        bounds_crs: str,
    ) -> xr.Dataset:
        """Subset dataset to spatial bounds."""
        west, south, east, north = bounds

        ds_crs_wkt = ds.spatial_ref.attrs.get("crs_wkt")
        if ds_crs_wkt is None:
            raise ValueError("Dataset missing CRS information")

        ds_crs = CRS.from_wkt(ds_crs_wkt)

        # Transform bounds to dataset CRS if needed
        if bounds_crs != ds_crs.to_string():
            left, bottom, right, top = rio_transform_bounds(
                CRS.from_string(bounds_crs),
                ds_crs,
                west,
                south,
                east,
                north,
            )
        else:
            left, bottom, right, top = west, south, east, north

        # Use latitude/longitude for subsetting
        ds_subset = ds.sel(
            longitude=slice(left, right),
            latitude=slice(top, bottom),
        )

        return ds_subset

    def _subset_height(self, ds: xr.Dataset, max_height: float) -> xr.Dataset:
        """Subset dataset by maximum height."""
        if "height" not in ds.dims:
            return ds

        return ds.where(ds["height"] <= max_height, drop=True)

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
        return f"TropoProduct(date={self.date.isoformat()}, model={self.model})"


# Functions


def interpolate_in_time(
    tropo_early: TropoProduct,
    tropo_late: TropoProduct,
    target_datetime: datetime,
    bounds: tuple[float, float, float, float] | None = None,
    max_height: float = 11e3,
    bounds_buffer: float = 0.2,
    output_path: Path | str | None = None,
) -> xr.DataArray:
    """Interpolate tropospheric delay between two products in time.

    Parameters
    ----------
    tropo_early : TropoProduct
        Earlier tropospheric product.
    tropo_late : TropoProduct
        Later tropospheric product.
    target_datetime : datetime
        Target datetime for interpolation.
    bounds : tuple[float, float, float, float] or None, optional
        Spatial bounds as (west, south, east, north). Default is None.
    max_height : float, optional
        Maximum height in meters. Default is 11000 m.
    bounds_buffer : float, optional
        Buffer to add to bounds in degrees. Default is 0.2.
    output_path : Path, str, or None, optional
        If provided, save result to NetCDF. Default is None.

    Returns
    -------
    xr.DataArray
        Interpolated tropospheric delay at target datetime.

    Examples
    --------
    >>> from datetime import datetime
    >>> from product import TropoProduct
    >>> from product._tropo import interpolate_in_time
    >>>
    >>> early = TropoProduct.from_path("tropo_00Z.nc")
    >>> late = TropoProduct.from_path("tropo_06Z.nc")
    >>> target = datetime(2022, 1, 11, 3, 0)
    >>>
    >>> # Basic interpolation
    >>> delay = interpolate_in_time(early, late, target)
    >>>
    >>> # With spatial subsetting
    >>> bounds = (-96, 29, -94, 31)
    >>> delay = interpolate_in_time(
    ...     early, late, target,
    ...     bounds=bounds,
    ...     max_height=11000,
    ...     bounds_buffer=0.2,
    ... )

    """
    if tropo_early.date > tropo_late.date:
        raise ValueError(
            f"Early product date ({tropo_early.date}) must be before "
            f"late product date ({tropo_late.date})"
        )

    if target_datetime < tropo_early.date or target_datetime > tropo_late.date:
        raise ValueError(
            f"Target datetime ({target_datetime}) must be between "
            f"early ({tropo_early.date}) and late ({tropo_late.date}) dates"
        )

    # Get total delay with consistent kwargs
    da_early = tropo_early.get_total_delay(
        bounds=bounds,
        max_height=max_height,
        bounds_buffer=bounds_buffer,
    )

    da_late = tropo_late.get_total_delay(
        bounds=bounds,
        max_height=max_height,
        bounds_buffer=bounds_buffer,
    )

    # Compute interpolation weight
    delta_total = (tropo_late.date - tropo_early.date).total_seconds()
    delta_target = (target_datetime - tropo_early.date).total_seconds()
    weight = delta_target / delta_total

    # Linear interpolation
    da_interp = (1 - weight) * da_early + weight * da_late

    # Add metadata
    da_interp.name = "zenith_total_delay"
    da_interp.attrs.update(
        {
            "interpolation_method": "linear",
            "early_product": tropo_early.filename,
            "late_product": tropo_late.filename,
            "early_date": tropo_early.date.isoformat(),
            "late_date": tropo_late.date.isoformat(),
            "target_date": target_datetime.isoformat(),
            "interpolation_weight": float(weight),
            "long_name": "Total zenith tropospheric delay",
            "units": "meters",
        }
    )

    # Save if requested
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        da_interp.to_netcdf(output_path, engine="h5netcdf")

    return da_interp


def interpolate_to_dem_surface(
    da_tropo_cube: xr.DataArray,
    dem: xr.DataArray,
    method: str = "linear",
    output_path: Path | str | None = None,
    output_format: str = "netcdf",
) -> xr.DataArray:
    """Interpolate 3D tropospheric delay to DEM surface heights.

    Parameters
    ----------
    da_tropo_cube : xr.DataArray
        3D tropospheric delay with dims (height, y, x).
        Assumed to be in EPSG:4326 (WGS84) if CRS not specified.
    dem : xr.DataArray
        DEM with surface heights in meters. Must have CRS information.
    method : str, optional
        Interpolation method ("linear" or "nearest"). Default is "linear".
    output_path : Path, str, or None, optional
        If provided, save result. Default is None.
    output_format : str, optional
        Output format ("netcdf" or "geotiff"). Default is "netcdf".

    Returns
    -------
    xr.DataArray
        2D tropospheric delay at DEM surface.

    Raises
    ------
    ValueError
        If DEM is missing CRS information.

    """
    # Ensure consistent coordinate naming
    if "latitude" in da_tropo_cube.dims:
        da_tropo_cube = da_tropo_cube.rename({"latitude": "y", "longitude": "x"})

    # Check DEM has CRS
    if not hasattr(dem, "rio") or dem.rio.crs is None:
        raise ValueError(
            "DEM is missing CRS information. "
            "Use dem.rio.write_crs() to set the CRS before calling this function."
        )

    dem_crs = dem.rio.crs

    # Write CRS to tropo if missing (assume EPSG:4326)
    if not hasattr(da_tropo_cube, "rio") or da_tropo_cube.rio.crs is None:
        da_tropo_cube = da_tropo_cube.rio.write_crs("EPSG:4326")

    # Reproject if different CRS
    if da_tropo_cube.rio.crs != dem_crs:
        td_utm = da_tropo_cube.rio.reproject(
            dem_crs,
            resampling=Resampling.cubic,
        )
    else:
        td_utm = da_tropo_cube

    # Find height dimension
    if "height" not in td_utm.dims:
        raise ValueError(
            f"No height dimension found. Available dims: {list(td_utm.dims)}"
        )

    # Build interpolator
    rgi = RegularGridInterpolator(
        (td_utm["height"].values, td_utm.y.values, td_utm.x.values),
        np.nan_to_num(td_utm.values),
        method=method,
        bounds_error=False,
        fill_value=np.nan,
    )

    # Create coordinate meshgrid
    yy, xx = np.meshgrid(dem.y.values, dem.x.values, indexing="ij")
    pts = np.column_stack([dem.values.ravel(), yy.ravel(), xx.ravel()])

    # Interpolate
    vals = rgi(pts)

    # Create output DataArray
    out = dem.copy()
    out.values[:] = vals.reshape(dem.shape).astype(np.float32)
    out.name = da_tropo_cube.name or "tropospheric_delay"

    # Update attributes
    out.attrs.update(
        {
            "interpolation_method": method,
            "interpolated_from": "3D tropospheric model",
            "units": "meters",
            "long_name": "Tropospheric delay at DEM surface",
        }
    )

    # Add time if present in input
    if "time" in td_utm.coords:
        out.attrs["time"] = str(td_utm.time.values)

    # Preserve any existing tropo metadata
    if "target_date" in da_tropo_cube.attrs:
        out.attrs["target_date"] = da_tropo_cube.attrs["target_date"]
    if "interpolation_weight" in da_tropo_cube.attrs:
        out.attrs["interpolation_weight"] = da_tropo_cube.attrs["interpolation_weight"]

    # Save if requested
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if output_format == "netcdf":
            out.to_netcdf(output_path, engine="h5netcdf")
        elif output_format == "geotiff":
            out.rio.to_raster(
                output_path,
                compress="deflate",
                tiled=True,
                dtype="float32",
            )
        else:
            raise ValueError(f"Unknown format: {output_format}")

    return out


def compute_los_correction(
    zenith_delay_2d: xr.DataArray,
    los_up: xr.DataArray,
    reference_correction: xr.DataArray | None = None,
    target_crs: str | None = None,
    output_path: Path | str | None = None,
    output_format: str = "geotiff",
) -> xr.DataArray:
    """Convert zenith delay to line-of-sight correction."""
    # Ensure same grid
    if los_up.shape != zenith_delay_2d.shape:
        if hasattr(los_up, "rio") and hasattr(zenith_delay_2d, "rio"):
            los_up = los_up.rio.reproject_match(zenith_delay_2d)
        else:
            raise ValueError(
                f"Shape mismatch: los_up {los_up.shape} vs "
                f"zenith_delay {zenith_delay_2d.shape}"
            )

    # Convert zenith to LOS
    # Note: -1 matches DISP convention (positive = apparent uplift)
    los_correction = -1 * (zenith_delay_2d / los_up)

    # Subtract reference if provided
    if reference_correction is not None:
        if reference_correction.shape != los_correction.shape:
            if hasattr(reference_correction, "rio"):
                reference_correction = reference_correction.rio.reproject_match(
                    los_correction
                )

        los_correction = (los_correction - reference_correction).astype(np.float32)

    # Reproject to target CRS if requested
    if target_crs is not None:
        if hasattr(los_correction, "rio"):
            los_correction = los_correction.rio.reproject(target_crs)
        else:
            raise ValueError("Cannot reproject: DataArray missing CRS information")

    # Add metadata
    los_correction.name = "los_correction"
    los_correction.attrs.update(
        {
            "units": "meters",
            "long_name": "Line-of-sight atmospheric correction",
            "line_of_sight_convention": (
                "Positive means decrease in delay (apparent uplift towards satellite)"
            ),
        }
    )

    if reference_correction is not None:
        los_correction.attrs["reference_subtracted"] = "yes"

    if target_crs is not None:
        los_correction.attrs["target_crs"] = target_crs

    # Save if requested
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if output_format == "netcdf":
            los_correction.to_netcdf(output_path, engine="h5netcdf")
        elif output_format == "geotiff":
            los_correction.rio.to_raster(
                output_path,
                compress="deflate",
                tiled=True,
                dtype="float32",
            )
        else:
            raise ValueError(f"Unknown format: {output_format}")

    return los_correction
