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

from ._disp import DispProduct


@dataclass
class CalProduct:
    """Calibrated OPERA DISP displacement product.

    Represents a calibration correction product that should be subtracted
    from OPERA DISP products. Main group contains calibration at full
    resolution. Optional model_3d group contains 3D displacement
    components at coarser resolution.

    Groups
    ------
    Main group:
        - calibration: Correction to subtract from DISP (full resolution)
        - calibration_std: Calibration uncertainty (full resolution)

    model_3d group (optional):
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
    ...     model_3d={"north_south": vel_ns, "east_west": vel_ew, "up_down": vel_ud},
    ... )

    >>> # Access main calibration
    >>> ds_main = cal.open_dataset()
    >>> calibration = ds_main["calibration"]

    >>> # Access 3D model (coarse resolution)
    >>> ds_model = cal.open_model_3d()
    >>> model_up = ds_model["up_down"]

    """

    path: Path
    frame_id: int
    primary_date: datetime
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
        r"(?P<primary>\d{8}T\d{6}Z)_"
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

        if self.secondary_date <= self.primary_date:
            raise ValueError(
                f"Secondary date ({self.secondary_date}) must be after "
                f"primary date ({self.primary_date})"
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
            primary_date=datetime.strptime(match.group("primary"), "%Y%m%dT%H%M%SZ"),
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
        disp_product: "DispProduct",
        output_dir: Path | str,
        sensor: str = "S1",
        calibration_std: xr.DataArray | None = None,
        model_3d: dict[str, xr.DataArray] | None = None,
        model_3d_std: dict[str, xr.DataArray] | None = None,
        spatial_ref: xr.DataArray | None = None,
        metadata: dict[str, str] | None = None,
        version: str = "1.0",
    ) -> "CalProduct":
        """Create calibration product with optional model_3d group.

        Parameters
        ----------
        calibration : xr.DataArray
            Calibration correction at full DISP resolution.
        disp_product : DispProduct
            Original DISP product (for metadata).
        output_dir : Path or str
            Output directory for NetCDF file. Filename auto-generated.
        sensor : str, optional
            Sensor type: "S1" or "NI". Default is "S1".
        calibration_std : xr.DataArray or None, optional
            Calibration uncertainty at full resolution. Default is None.
        model_3d : dict[str, xr.DataArray] or None, optional
            3D displacement components (coarse resolution) with keys:
            "north_south", "east_west", "up_down". Default is None.
        model_3d_std : dict[str, xr.DataArray] or None, optional
            3D displacement uncertainties (coarse resolution). Default is None.
        spatial_ref : xr.DataArray or None, optional
            Spatial reference data variable from input DISP product. Default is None.
        metadata : dict[str, str] or None, optional
            Additional metadata (e.g., GNSS reference epoch). Default is None.
        version : str, optional
            Product version. Default is "1.0".

        Returns
        -------
        CalProduct
            Created calibration product.

        Examples
        --------
        >>> from product import DispProduct, CalProduct
        >>>
        >>> disp = DispProduct.from_path("input/disp.nc")
        >>> ds = disp.open_dataset()
        >>>
        >>> cal = CalProduct.create(
        ...     calibration=cal_full,
        ...     disp_product=disp,
        ...     output_dir="output/",
        ...     calibration_std=cal_std,
        ...     spatial_ref=ds["spatial_ref"],
        ...     metadata={"gnss_reference_epoch": "2020-01-01T00:00:00Z"},
        ... )
        >>> print(cal.filename)
        OPERA_L4_CAL-DISP-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0_20260104T123456Z.nc

        """
        import rioxarray  # noqa: F401

        if sensor not in {"S1", "NI"}:
            raise ValueError(f"Invalid sensor: {sensor}. Must be 'S1' or 'NI'")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate OPERA-compliant filename
        production_date = datetime.utcnow()
        filename = (
            f"OPERA_L4_CAL-DISP-{sensor}_"
            f"{disp_product.mode}_"
            f"F{disp_product.frame_id:05d}_"
            f"{disp_product.polarization}_"
            f"{disp_product.primary_date:%Y%m%dT%H%M%S}Z_"
            f"{disp_product.secondary_date:%Y%m%dT%H%M%S}Z_"
            f"v{version}_"
            f"{production_date:%Y%m%dT%H%M%S}Z.nc"
        )
        output_file = output_dir / filename

        # Build main group dataset (full resolution)
        data_vars = {"calibration": calibration}

        if calibration_std is not None:
            data_vars["calibration_std"] = calibration_std

        # Add spatial_ref as data variable if provided
        if spatial_ref is not None:
            data_vars["spatial_ref"] = spatial_ref

        # Create main dataset
        ds = xr.Dataset(data_vars)

        # Write CRS using rioxarray if spatial_ref contains crs_wkt
        if spatial_ref is not None:
            crs_wkt = spatial_ref.attrs.get("crs_wkt")
            if crs_wkt:
                ds = ds.rio.write_crs(crs_wkt)

        # Add global attributes
        ds.attrs.update(
            {
                "product_type": f"OPERA_L4_CAL-DISP-{sensor}",
                "sensor": sensor,
                "frame_id": disp_product.frame_id,
                "mode": disp_product.mode,
                "polarization": disp_product.polarization,
                "primary_datetime": disp_product.primary_date.isoformat(),
                "secondary_datetime": disp_product.secondary_date.isoformat(),
                "production_datetime": production_date.isoformat(),
                "product_version": version,
                "description": (
                    f"Calibration correction for {sensor} InSAR displacement (subtract"
                    " from DISP)"
                ),
                "source_product": disp_product.filename,
                "usage": (
                    "Subtract calibration layer from DISP displacement to obtain"
                    " calibrated displacement"
                ),
            }
        )

        # Add custom metadata
        if metadata:
            ds.attrs.update(metadata)

        # Save main group
        ds.to_netcdf(output_file, engine="h5netcdf")

        # Create model_3d group if 3D components provided (coarse resolution)
        if model_3d or model_3d_std:
            model_data_vars = {}

            # Add 3D displacement components
            if model_3d:
                for comp in ["north_south", "east_west", "up_down"]:
                    if comp in model_3d:
                        model_data_vars[comp] = model_3d[comp]

            # Add 3D displacement uncertainties
            if model_3d_std:
                for comp in ["north_south_std", "east_west_std", "up_down_std"]:
                    if comp in model_3d_std:
                        model_data_vars[comp] = model_3d_std[comp]

            # Add spatial_ref to model_3d group as well
            if spatial_ref is not None:
                model_data_vars["spatial_ref"] = spatial_ref

            if model_data_vars:
                ds_model = xr.Dataset(model_data_vars)

                # Write CRS to model_3d group if available
                if spatial_ref is not None:
                    crs_wkt = spatial_ref.attrs.get("crs_wkt")
                    if crs_wkt:
                        ds_model = ds_model.rio.write_crs(crs_wkt)

                # Add model group attributes
                ds_model.attrs.update(
                    {
                        "description": (
                            "3D displacement model at coarse resolution (e.g., from"
                            " GNSS grid or deformation model)"
                        ),
                        "units": "meters",
                        "reference_frame": "ENU (East-North-Up)",
                    }
                )

                # Append to existing file as model_3d group
                ds_model.to_netcdf(
                    output_file,
                    mode="a",
                    group="model_3d",
                    engine="h5netcdf",
                )

        return cls(
            path=output_file,
            frame_id=disp_product.frame_id,
            primary_date=disp_product.primary_date,
            secondary_date=disp_product.secondary_date,
            polarization=disp_product.polarization,
            sensor=sensor,
            version=version,
            production_date=production_date,
            mode=disp_product.mode,
        )

    def open_dataset(self, group: str | None = None) -> xr.Dataset:
        """Open calibration dataset.

        Parameters
        ----------
        group : str or None, optional
            Group to open: None for main, "model_3d" for 3D model.
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

        >>> # Open model_3d group (coarse resolution)
        >>> ds_model = cal.open_dataset(group="model_3d")
        >>> model_up = ds_model["up_down"]

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        if group == "model_3d":
            return xr.open_dataset(self.path, group="model_3d", engine="h5netcdf")

        return xr.open_dataset(self.path, engine="h5netcdf")

    def open_model_3d(self) -> xr.Dataset:
        """Open model_3d group dataset.

        Returns
        -------
        xr.Dataset
            Dataset containing 3D displacement model at coarse resolution.

        Raises
        ------
        FileNotFoundError
            If product file does not exist.
        ValueError
            If model_3d group does not exist.

        Examples
        --------
        >>> ds_model = cal.open_model_3d()
        >>> disp_ns = ds_model["north_south"]
        >>> disp_ew = ds_model["east_west"]
        >>> disp_up = ds_model["up_down"]

        """
        if not self.path.exists():
            raise FileNotFoundError(f"Product file not found: {self.path}")

        try:
            return xr.open_dataset(self.path, group="model_3d", engine="h5netcdf")
        except (OSError, ValueError) as e:
            raise ValueError(
                f"model_3d group not found in {self.filename}. "
                "Product may not contain 3D displacement model."
            ) from e

    def has_model_3d(self) -> bool:
        """Check if product contains model_3d group.

        Returns
        -------
        bool
            True if model_3d group exists.

        """
        try:
            self.open_model_3d()
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
        group : str or None, optional
            Group containing layer. Default is None (main group).
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
        >>> # Export main calibration
        >>> cal.to_geotiff("calibration", "calibration.tif")

        >>> # Export 3D model component
        >>> cal.to_geotiff("up_down", "model_up.tif", group="model_3d")

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
        data = da.values

        # Extract spatial information
        if "spatial_ref" in ds:
            transform = self._get_transform(ds)
            crs = self._get_crs(ds)
        else:
            transform = Affine.translation(
                float(ds.x.values[0]),
                float(ds.y.values[0]),
            ) * Affine.scale(
                float(ds.x.values[1] - ds.x.values[0]),
                float(ds.y.values[1] - ds.y.values[0]),
            )
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
                primary_date=self.primary_date.isoformat(),
                secondary_date=self.secondary_date.isoformat(),
                layer=layer,
                group=group if group else "main",
            )

        return output_path

    def get_calibration_summary(self) -> dict[str, dict[str, float]]:
        """Get summary statistics of all layers.

        Returns
        -------
        dict[str, dict[str, float]]
            Statistics for main and model_3d groups.

        Examples
        --------
        >>> summary = cal.get_calibration_summary()
        >>> summary["main"]["calibration"]
        {'mean': 0.023, 'std': 0.015, 'min': -0.05, 'max': 0.08}
        >>> summary["model_3d"]["up_down"]
        {'mean': 0.001, 'std': 0.003, 'min': -0.01, 'max': 0.02}

        """
        summary: dict = {"main": {}}

        # Main group
        ds = self.open_dataset()
        for var in ds.data_vars:
            if var == "spatial_ref":
                continue

            data = ds[var].values
            valid_data = data[~np.isnan(data)]

            if len(valid_data) > 0:
                summary["main"][var] = {
                    "mean": float(np.mean(valid_data)),
                    "std": float(np.std(valid_data)),
                    "min": float(np.min(valid_data)),
                    "max": float(np.max(valid_data)),
                }

        # model_3d group if exists
        if self.has_model_3d():
            summary["model_3d"] = {}
            ds_model = self.open_model_3d()

            for var in ds_model.data_vars:
                if var == "spatial_ref":
                    continue

                data = ds_model[var].values
                valid_data = data[~np.isnan(data)]

                if len(valid_data) > 0:
                    summary["model_3d"][var] = {
                        "mean": float(np.mean(valid_data)),
                        "std": float(np.std(valid_data)),
                        "min": float(np.min(valid_data)),
                        "max": float(np.max(valid_data)),
                    }

        return summary

    def _get_transform(self, ds: xr.Dataset) -> Affine:
        """Extract affine transform from dataset.

        Parameters
        ----------
        ds : xr.Dataset
            Dataset with spatial_ref containing GeoTransform.

        Returns
        -------
        Affine
            Affine transformation.

        Raises
        ------
        ValueError
            If GeoTransform not found.

        """
        gt = ds.spatial_ref.attrs.get("GeoTransform")
        if gt is None:
            raise ValueError("No GeoTransform found in spatial_ref")

        vals = [float(x) for x in gt.split()]
        return Affine(vals[1], vals[2], vals[0], vals[4], vals[5], vals[3])

    def _get_crs(self, ds: xr.Dataset) -> str:
        """Extract CRS from dataset.

        Parameters
        ----------
        ds : xr.Dataset
            Dataset with spatial_ref containing crs_wkt.

        Returns
        -------
        str
            CRS as WKT string.

        Raises
        ------
        ValueError
            If crs_wkt not found.

        """
        crs_wkt = ds.spatial_ref.attrs.get("crs_wkt")
        if crs_wkt is None:
            raise ValueError("No crs_wkt found in spatial_ref")
        return crs_wkt

    @property
    def baseline_days(self) -> int:
        """Temporal baseline in days."""
        return (self.secondary_date - self.primary_date).days

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
        model_str = "+model_3d" if self.exists and self.has_model_3d() else ""
        return (
            f"CalProduct(sensor={self.sensor}, frame={self.frame_id}, "
            f"{self.primary_date.date()} → {self.secondary_date.date()}, "
            f"{self.polarization}{model_str})"
        )
