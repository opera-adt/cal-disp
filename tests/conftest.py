"""Pytest configuration and shared fixtures for cal-disp."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Iterator

import numpy as np
import pytest
import xarray as xr
from numpy.typing import NDArray

# Basic directory fixtures


@pytest.fixture
def temp_work_dir(tmp_path: Path) -> Path:
    """Temporary working directory."""
    work_dir = tmp_path / "work"
    work_dir.mkdir()
    return work_dir


@pytest.fixture
def temp_output_dir(tmp_path: Path) -> Path:
    """Temporary output directory."""
    output_dir = tmp_path / "outputs"
    output_dir.mkdir()
    return output_dir


@pytest.fixture
def sample_frame_id() -> int:
    """Sample OPERA frame ID."""
    return 8882


@pytest.fixture
def sample_grid_version() -> str:
    """Sample UNR gridded data version."""
    return "0.2"


@pytest.fixture
def sample_grid_type() -> str:
    """Sample UNR gridded data type."""
    return "constant"


# DISP product fixtures


@pytest.fixture
def sample_disp_product(tmp_path: Path) -> Path:
    """Create a mock DISP-S1 product file."""
    ny, nx = 100, 100

    ds = xr.Dataset(
        {
            "displacement": (["y", "x"], np.random.randn(ny, nx)),
            "temporal_coherence": (["y", "x"], np.random.uniform(0.3, 0.9, (ny, nx))),
            "conncomp": (["y", "x"], np.ones((ny, nx), dtype=np.uint16)),
            "latitude": (
                ["y", "x"],
                np.linspace(34.0, 35.0, ny)[:, None] * np.ones((ny, nx)),
            ),
            "longitude": (
                ["y", "x"],
                np.linspace(-118.0, -117.0, nx)[None, :] * np.ones((ny, nx)),
            ),
        },
        coords={
            "y": np.arange(ny),
            "x": np.arange(nx),
            "time": [datetime(2022, 7, 22, 0, 26, 57)],
        },
        attrs={
            "description": "OPERA DISP-S1 Product",
            "institution": "NASA JPL",
            "mission_name": "OPERA",
            "reference_datetime": "20220111T002651Z",
            "secondary_datetime": "20220722T002657Z",
        },
    )

    filename = (
        "OPERA_L3_DISP-S1_IW_F08882_VV_"
        "20220111T002651Z_20220722T002657Z_"
        "v1.0_20230101T000000Z.nc"
    )
    filepath = tmp_path / filename
    ds.to_netcdf(filepath, engine="h5netcdf")
    return filepath


@pytest.fixture
def sample_disp_product_with_corrections(tmp_path: Path) -> Path:
    """Create mock DISP product with corrections group."""
    ny, nx = 50, 50

    # Main group with spatial_ref
    spatial_ref = xr.DataArray(
        0,
        attrs={
            "crs_wkt": (
                'PROJCS["WGS 84 / UTM zone 15N",'
                'GEOGCS["WGS 84",DATUM["WGS_1984",'
                'SPHEROID["WGS 84",6378137,298.257223563]],'
                'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
                'PROJECTION["Transverse_Mercator"],'
                'PARAMETER["latitude_of_origin",0],'
                'PARAMETER["central_meridian",-93],'
                'PARAMETER["scale_factor",0.9996],'
                'PARAMETER["false_easting",500000],'
                'PARAMETER["false_northing",0],UNIT["metre",1]]'
            ),
            "GeoTransform": "100000.0 30.0 0.0 3500000.0 0.0 -30.0",
        },
    )

    ds = xr.Dataset(
        {
            "displacement": (["y", "x"], np.random.randn(ny, nx) * 0.1),
            "temporal_coherence": (["y", "x"], np.random.rand(ny, nx)),
            "spatial_ref": spatial_ref,
        },
        coords={
            "y": np.arange(ny) * 30.0 + 3500000.0,
            "x": np.arange(nx) * 30.0 + 100000.0,
            "time": [datetime(2022, 7, 22, 0, 26, 57)],
        },
    )

    # Corrections group
    ref_point = xr.DataArray(
        0,
        attrs={
            "rows": 25,
            "cols": 25,
            "latitudes": 34.5,
            "longitudes": -118.0,
        },
    )

    ds_corr = xr.Dataset(
        {
            "ionospheric_delay": (["y", "x"], np.random.randn(ny, nx) * 0.01),
            "solid_earth_tide": (["y", "x"], np.random.randn(ny, nx) * 0.001),
            "reference_point": ref_point,
        },
        coords={
            "y": np.arange(ny) * 30.0 + 3500000.0,
            "x": np.arange(nx) * 30.0 + 100000.0,
        },
    )

    filename = (
        "OPERA_L3_DISP-S1_IW_F08882_VV_"
        "20220111T002651Z_20220722T002657Z_"
        "v1.0_20230101T000000Z.nc"
    )
    filepath = tmp_path / filename

    ds.to_netcdf(filepath, engine="h5netcdf")
    ds_corr.to_netcdf(filepath, group="corrections", mode="a", engine="h5netcdf")

    return filepath


# Static layer fixtures


@pytest.fixture
def sample_static_layers(tmp_path: Path) -> tuple[Path, Path]:
    """Create mock LOS and DEM static layer files."""
    ny, nx = 100, 100

    # LOS file
    los_ds = xr.Dataset(
        {
            "los_east": (["y", "x"], np.random.randn(ny, nx) * 0.1),
            "los_north": (["y", "x"], np.random.randn(ny, nx) * 0.1),
            "los_up": (["y", "x"], np.ones((ny, nx)) * 0.9),
        },
        coords={"y": np.arange(ny), "x": np.arange(nx)},
    )
    los_file = tmp_path / "los.h5"
    los_ds.to_netcdf(los_file, engine="h5netcdf")

    # DEM file
    dem_ds = xr.Dataset(
        {
            "elevation": (["y", "x"], np.random.uniform(0, 1000, (ny, nx))),
        },
        coords={"y": np.arange(ny), "x": np.arange(nx)},
    )
    dem_file = tmp_path / "dem.h5"
    dem_ds.to_netcdf(dem_file, engine="h5netcdf")

    return los_file, dem_file


@pytest.fixture
def sample_static_dem(tmp_path: Path) -> Path:
    """Create mock DEM static layer GeoTIFF."""
    try:
        import rasterio
        from rasterio.crs import CRS
        from rasterio.transform import from_bounds

        dem_data = np.random.uniform(0, 2000, (100, 100)).astype(np.float32)
        dem_file = tmp_path / "OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_dem.tif"

        transform = from_bounds(-118, 34, -117, 35, 100, 100)
        crs = CRS.from_epsg(4326)

        with rasterio.open(
            dem_file,
            "w",
            driver="GTiff",
            height=100,
            width=100,
            count=1,
            dtype=np.float32,
            crs=crs,
            transform=transform,
            nodata=-9999,
        ) as dst:
            dst.write(dem_data, 1)

        return dem_file
    except ImportError:
        pytest.skip("rasterio not available")
        return tmp_path / "dummy.tif"


@pytest.fixture
def sample_static_los(tmp_path: Path) -> Path:
    """Create mock LOS static layer GeoTIFF with 3 bands."""
    try:
        import rasterio
        from rasterio.crs import CRS
        from rasterio.transform import from_bounds

        # Create LOS unit vectors (east, north, up)
        los_east = np.random.randn(100, 100).astype(np.float32) * 0.1
        los_north = np.random.randn(100, 100).astype(np.float32) * 0.1
        los_up = np.ones((100, 100), dtype=np.float32) * 0.9

        los_file = (
            tmp_path
            / "OPERA_L3_DISP-S1-STATIC_F08882_20140403_S1A_v1.0_line_of_sight_enu.tif"
        )

        transform = from_bounds(-118, 34, -117, 35, 100, 100)
        crs = CRS.from_epsg(4326)

        with rasterio.open(
            los_file,
            "w",
            driver="GTiff",
            height=100,
            width=100,
            count=3,
            dtype=np.float32,
            crs=crs,
            transform=transform,
            nodata=-9999,
        ) as dst:
            dst.write(los_east, 1)
            dst.write(los_north, 2)
            dst.write(los_up, 3)
            dst.set_band_description(1, "LOS East")
            dst.set_band_description(2, "LOS North")
            dst.set_band_description(3, "LOS Up")

        return los_file
    except ImportError:
        pytest.skip("rasterio not available")
        return tmp_path / "dummy.tif"


# Tropospheric product fixtures


@pytest.fixture
def sample_tropo_product(tmp_path: Path) -> Path:
    """Create mock TROPO-ZENITH product."""
    nh, nlat, nlon = 20, 50, 50

    # Create 3D grids
    height = np.linspace(0, 15000, nh)
    lat = np.linspace(34, 35, nlat)
    lon = np.linspace(-118, -117, nlon)

    # Random tropospheric delays
    wet = np.random.rand(nh, nlat, nlon) * 0.3
    hydro = np.random.rand(nh, nlat, nlon) * 2.3

    spatial_ref = xr.DataArray(
        0,
        attrs={
            "crs_wkt": (
                'GEOGCS["WGS 84",DATUM["WGS_1984",'
                'SPHEROID["WGS 84",6378137,298.257223563]],'
                'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]'
            ),
        },
    )

    ds = xr.Dataset(
        {
            "wet_delay": (["height", "latitude", "longitude"], wet),
            "hydrostatic_delay": (["height", "latitude", "longitude"], hydro),
            "spatial_ref": spatial_ref,
        },
        coords={
            "height": height,
            "latitude": lat,
            "longitude": lon,
            "time": [datetime(2022, 1, 11, 0, 0)],
        },
    )

    filename = "OPERA_L4_TROPO-ZENITH_20220111T000000Z_20220111T120000Z_HRRR_v1.0.nc"
    filepath = tmp_path / filename
    ds.to_netcdf(filepath, engine="h5netcdf")

    return filepath


# UNR GNSS fixtures


@pytest.fixture
def sample_unr_grid_latlon(tmp_path: Path) -> Path:
    """Create sample UNR grid lookup table."""
    content = """# UNR Grid Lookup Table v0.2
# lat lon grid_id
34.0 -118.0 1
34.1 -118.0 2
34.0 -117.9 3
34.1 -117.9 4
"""
    grid_file = tmp_path / "grid_latlon_lookup_v0.2.txt"
    grid_file.write_text(content)
    return grid_file


@pytest.fixture
def sample_unr_timeseries_dir(tmp_path: Path) -> Path:
    """Create directory with sample UNR .tenv8 files."""
    ts_dir = tmp_path / "unr_ts"
    ts_dir.mkdir()

    # Create a few mock .tenv8 files
    for i in range(1, 4):
        ts_file = ts_dir / f"{i:06d}_TEST.tenv8"
        content = f"""# UNR Timeseries Station {i}
# YYMMMDD       __east(mm)        __north(mm)       ____up(mm)
2022.0000    0.0    0.0    0.0    1.0    1.0    1.0  0
2022.0833    1.2    0.5    2.1    1.0    1.0    1.0  0
2022.1667    2.4    1.0    4.2    1.0    1.0    1.0  0
"""
        ts_file.write_text(content)

    return ts_dir


@pytest.fixture
def sample_unr_data(tmp_path: Path) -> tuple[Path, Path]:
    """Create sample UNR lookup table and tenv8 files."""
    # Lookup table - format: grid_point lon lat (not lat lon!)
    lookup_content = """000001 -118.0 34.0
000002 -118.0 34.1
000003 -117.9 34.0
000004 -117.9 34.1
"""
    lookup_file = tmp_path / "grid_latlon_lookup_v0.2.txt"
    lookup_file.write_text(lookup_content)

    # tenv8 files directory
    tenv8_dir = tmp_path / "tenv8"
    tenv8_dir.mkdir()

    # Create sample tenv8 files
    for grid_id in [1, 2, 3, 4]:
        tenv8_file = tenv8_dir / f"{grid_id:06d}_TEST.tenv8"
        content = """2022.0000    0.0    0.0    0.0    1.0    1.0    1.0  0
2022.0833    1.2    0.5    2.1    1.0    1.0    1.0  0
2022.1667    2.4    1.0    4.2    1.0    1.0    1.0  0
"""
        tenv8_file.write_text(content)

    return lookup_file, tenv8_dir


# Algorithm and configuration fixtures


@pytest.fixture
def sample_algorithm_params(tmp_path: Path) -> Path:
    """Create sample algorithm parameters YAML file."""
    params_content = """
algorithm_parameters:
  grid_search:
    enabled: true
    min_stations: 3
    max_distance_km: 100

  robust_estimation:
    method: "huber"
    outlier_threshold: 3.0

  quality_control:
    min_coherence: 0.3
    max_phase_std: 2.0
"""
    params_file = tmp_path / "algorithm_params.yaml"
    params_file.write_text(params_content)
    return params_file


@pytest.fixture
def sample_mask_file(tmp_path: Path) -> Path:
    """Create sample byte mask file."""
    try:
        import rasterio
        from rasterio.transform import from_bounds

        mask_data = np.random.randint(0, 2, (100, 100), dtype=np.uint8)
        mask_file = tmp_path / "mask.tif"
        transform = from_bounds(-118, 34, -117, 35, 100, 100)

        with rasterio.open(
            mask_file,
            "w",
            driver="GTiff",
            height=100,
            width=100,
            count=1,
            dtype=np.uint8,
            transform=transform,
        ) as dst:
            dst.write(mask_data, 1)

        return mask_file
    except ImportError:
        # If rasterio not available, create dummy file
        mask_file = tmp_path / "mask.tif"
        mask_file.touch()
        return mask_file


# Mock API responses


@pytest.fixture
def mock_cmr_response() -> list[dict]:
    """Mock CMR API response for DISP products."""
    return [
        {
            "producer_granule_id": (
                "OPERA_L3_DISP-S1_IW_F08882_VV_20220111T002651Z_"
                "20220123T002651Z_v1.0.nc"
            ),
            "time_start": "2022-01-11T00:26:51.000Z",
            "time_end": "2022-01-23T00:26:51.000Z",
            "granule_size": "123456789",
            "links": [
                {
                    "rel": "http://esipfed.org/ns/fedsearch/1.1/data#",
                    "href": "https://fake-url.earthdata.nasa.gov/product1.nc",
                }
            ],
        },
        {
            "producer_granule_id": (
                "OPERA_L3_DISP-S1_IW_F08882_VV_20220123T002651Z_"
                "20220204T002651Z_v1.0.nc"
            ),
            "time_start": "2022-01-23T00:26:51.000Z",
            "time_end": "2022-02-04T00:26:51.000Z",
            "granule_size": "123456789",
            "links": [
                {
                    "rel": "http://esipfed.org/ns/fedsearch/1.1/data#",
                    "href": "https://fake-url.earthdata.nasa.gov/product2.nc",
                }
            ],
        },
    ]


# Synthetic data fixtures


@pytest.fixture
def sample_date_range() -> tuple[datetime, datetime]:
    """Sample date range for testing."""
    start = datetime(2023, 1, 1)
    end = datetime(2023, 1, 31)
    return start, end


@pytest.fixture
def sample_phase_data() -> NDArray[np.float64]:
    """Generate synthetic unwrapped phase data."""
    rng = np.random.default_rng(42)
    ny, nx = 100, 100

    # Create some spatial correlation
    x = np.linspace(0, 2 * np.pi, nx)
    y = np.linspace(0, 2 * np.pi, ny)
    X, Y = np.meshgrid(x, y)

    phase = 5.0 * np.sin(X) * np.cos(Y) + rng.normal(0, 0.5, (ny, nx))

    return phase


@pytest.fixture
def sample_coherence() -> NDArray[np.float64]:
    """Generate synthetic coherence data."""
    rng = np.random.default_rng(42)
    ny, nx = 100, 100

    # Higher coherence in center, lower at edges
    y = np.linspace(-1, 1, ny)
    x = np.linspace(-1, 1, nx)
    X, Y = np.meshgrid(x, y)

    # Gaussian-like coherence pattern
    coherence = 0.9 * np.exp(-(X**2 + Y**2) / 0.5)
    coherence += rng.uniform(0, 0.1, (ny, nx))
    coherence = np.clip(coherence, 0, 1)

    return coherence


@pytest.fixture
def sample_coordinates() -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Generate sample lat/lon coordinates."""
    ny, nx = 100, 100

    # California-ish coordinates
    lat = np.linspace(34.0, 35.0, ny)
    lon = np.linspace(-120.0, -119.0, nx)

    lon_grid, lat_grid = np.meshgrid(lon, lat)

    return lat_grid, lon_grid


@pytest.fixture
def sample_disp_dataset(
    sample_phase_data: NDArray[np.float64],
    sample_coherence: NDArray[np.float64],
    sample_coordinates: tuple[NDArray[np.float64], NDArray[np.float64]],
) -> xr.Dataset:
    """Create a synthetic DISP-S1 dataset."""
    lat, lon = sample_coordinates
    ny, nx = lat.shape

    y = np.arange(ny)
    x = np.arange(nx)

    ds = xr.Dataset(
        {
            "unwrapped_phase": (["y", "x"], sample_phase_data),
            "coherence": (["y", "x"], sample_coherence),
            "latitude": (["y", "x"], lat),
            "longitude": (["y", "x"], lon),
            "solid_earth_tide": (["y", "x"], np.zeros((ny, nx))),
            "ionosphere": (["y", "x"], np.zeros((ny, nx))),
            "troposphere": (["y", "x"], np.zeros((ny, nx))),
        },
        coords={
            "y": y,
            "x": x,
        },
        attrs={
            "description": "OPERA DISP-S1 Product",
            "institution": "NASA JPL",
            "mission_name": "OPERA",
            "reference_datetime": "2023-01-01T00:00:00",
            "secondary_datetime": "2023-01-13T00:00:00",
        },
    )

    return ds


@pytest.fixture
def sample_disp_file(
    tmp_path: Path,
    sample_disp_dataset: xr.Dataset,
) -> Path:
    """Create a temporary NetCDF file with DISP data."""
    filepath = tmp_path / "test_disp_product.nc"
    sample_disp_dataset.to_netcdf(filepath, engine="h5netcdf")
    return filepath


# Pytest configuration


@pytest.fixture
def cli_runner():
    """Click CLI test runner."""
    from click.testing import CliRunner

    return CliRunner()


@pytest.fixture(autouse=True)
def reset_random_state() -> Iterator[None]:
    """Reset numpy random state before each test."""
    np.random.seed(42)
    yield


def pytest_configure(config):
    """Add custom markers."""
    config.addinivalue_line(
        "markers",
        "slow: marks tests as slow (deselect with '-m \"not slow\"')",
    )
    config.addinivalue_line(
        "markers",
        "integration: integration tests requiring external data or services",
    )
    config.addinivalue_line(
        "markers",
        "requires_earthdata: tests requiring Earthdata credentials",
    )
    config.addinivalue_line(
        "markers",
        "requires_network: tests requiring network access",
    )
