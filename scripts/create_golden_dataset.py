"""Generate a synthetic golden dataset for cal-disp validation.

Creates small, fully deterministic input files that mirror the schema of real
OPERA DISP-S1 / F36540 products, runs the complete calibration workflow, and
saves both inputs and expected output to an external directory.

The generated output can be used as a reference with ``cal-disp validate``
to check that future runs on the same inputs produce identical results:

    python scripts/create_golden_dataset.py --output-dir /path/to/test_data
    cal-disp run new_runconfig.yaml
    cal-disp validate /path/to/test_data/golden_output/<ref>.nc <new_output.nc>

For automated workflow testing (no golden data required) use:

    pytest -m integration

Outputs (written under <output-dir>/golden_datasets/)
------------------------------------------------------
input_data/
    disp/                              -- 200 × 200 synthetic DISP-S1 NetCDF
                                          (30 m, UTM 11N)
    gnss/                              -- UNR lookup table + 4 × .tenv8 files
    static_input/
        OPERA_L3_DISP-S1-STATIC_F*_line_of_sight_enu.tif
        OPERA_L3_DISP-S1-STATIC_F*_dem.tif
    tropo/
        OPERA_L4_TROPO-ZENITH_<ref_date>_..._HRES_v1.0.nc
        OPERA_L4_TROPO-ZENITH_<sec_date>_..._HRES_v1.0.nc

configs/
    algorithm_parameters.yaml

golden_output/
    OPERA_L4_DISP-CAL-S1_IW_F36540_VV_*.nc  -- known-good reference CalProduct

output/            -- empty; reserved for new runs during validation

Re-running the script regenerates all files deterministically.
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

REPO_ROOT = Path(__file__).resolve().parent.parent
_ENV_VAR = "CAL_DISP_TEST_DATA"


def _resolve_output_dir(cli_arg: str | None) -> Path:
    """Return the output root, preferring the CLI arg over the env var."""
    if cli_arg:
        return Path(cli_arg).resolve()
    env = os.environ.get(_ENV_VAR)
    if env:
        return Path(env).resolve()
    sys.exit(
        f"ERROR: Provide --output-dir or set the ${_ENV_VAR} environment variable."
    )


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--output-dir",
        metavar="DIR",
        default=None,
        help=f"Root directory for the golden dataset (default: ${_ENV_VAR} env var).",
    )
    return p.parse_args()


# Resolved at call-time in main(); module-level vars set there.
INPUT_DATA_DIR: Path
CONFIGS_DIR: Path
GOLDEN_OUTPUT_DIR: Path

# Grid constants — match real frame-36540 resolution and CRS
NY, NX = 200, 200
SPACING = 30.0  # meters (same as real data)
EPSG = 32611  # WGS 84 / UTM zone 11N

# Top-left pixel *centre* coords (UTM 11N meters)
# Chosen so that synthetic GNSS stations sit inside the extent.
X0 = 405000.0
Y0 = 3778000.0

# Acquisition dates
REF_DATETIME = datetime(2016, 7, 24, 1, 58, 9, tzinfo=timezone.utc)
SEC_DATETIME = datetime(2016, 8, 5, 1, 58, 9, tzinfo=timezone.utc)
FRAME_ID = 36540
DISP_FILENAME = (
    f"OPERA_L3_DISP-S1_IW_F{FRAME_ID:05d}_VV_"
    f"{REF_DATETIME:%Y%m%dT%H%M%S}Z_"
    f"{SEC_DATETIME:%Y%m%dT%H%M%S}Z_"
    "v1.0_20250101T000000Z.nc"
)
LOS_FILENAME = (
    f"OPERA_L3_DISP-S1-STATIC_F{FRAME_ID:05d}_20140403_S1A_v1.0_line_of_sight_enu.tif"
)
DEM_FILENAME = f"OPERA_L3_DISP-S1-STATIC_F{FRAME_ID:05d}_20140403_S1A_v1.0_dem.tif"

# Tropo filenames
_TROPO_PROD = "20250101T000000Z"
REF_TROPO_FILENAME = (
    f"OPERA_L4_TROPO-ZENITH_{REF_DATETIME:%Y%m%dT%H%M%S}Z_{_TROPO_PROD}_HRES_v1.0.nc"
)
SEC_TROPO_FILENAME = (
    f"OPERA_L4_TROPO-ZENITH_{SEC_DATETIME:%Y%m%dT%H%M%S}Z_{_TROPO_PROD}_HRES_v1.0.nc"
)

RNG = np.random.default_rng(0)


# Helper: CRS WKT for EPSG:32611
def _utm11n_crs_wkt() -> str:
    try:
        from pyproj import CRS

        return CRS.from_epsg(EPSG).to_wkt()
    except ImportError:
        # Fallback: hard-coded compact WKT matching real product
        return (
            'PROJCS["WGS 84 / UTM zone 11N",'
            'GEOGCS["WGS 84",DATUM["WGS_1984",'
            'SPHEROID["WGS 84",6378137,298.257223563]],'
            'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
            'PROJECTION["Transverse_Mercator"],'
            'PARAMETER["latitude_of_origin",0],'
            'PARAMETER["central_meridian",-117],'
            'PARAMETER["scale_factor",0.9996],'
            'PARAMETER["false_easting",500000],'
            'PARAMETER["false_northing",0],'
            'UNIT["metre",1],AUTHORITY["EPSG","32611"]]'
        )


def _utm_to_wgs84(easting: float, northing: float) -> tuple[float, float]:
    """Return (lon, lat) in WGS84 for a UTM 11N (east, north) point."""
    try:
        from pyproj import Transformer

        t = Transformer.from_crs(EPSG, 4326, always_xy=True)
        return t.transform(easting, northing)
    except ImportError:
        # Approximate conversion (good to ~0.01° for this area)
        false_e, cm = 500000.0, -117.0
        k0, a = 0.9996, 6378137.0
        n = (northing) / (k0 * a)
        lon = cm + (easting - false_e) / (k0 * a * np.cos(np.radians(34.0))) * (
            180 / np.pi
        )
        lat = n * (180 / np.pi)
        return lon, lat


# Golden DISP NetCDF
def create_disp(out_dir: Path) -> Path:
    """Create a synthetic DISP-S1 NetCDF product and return its path."""
    out_dir.mkdir(parents=True, exist_ok=True)

    x = np.arange(NX, dtype=np.float64) * SPACING + X0
    y = np.arange(NY, dtype=np.float64) * (-SPACING) + Y0

    # Planar ramp (long-wavelength signal typical of troposphere / orbital error)
    xx, yy = np.meshgrid(
        (x - x.mean()) / (NX * SPACING),
        (y - y.mean()) / (NY * SPACING),
    )
    ramp = 0.005 * xx + 0.003 * yy  # meters

    # Gaussian blob simulating localised deformation
    r2 = xx**2 + yy**2
    blob = 0.012 * np.exp(-r2 / 0.1)

    displacement = (ramp + blob).astype(np.float32)

    # Temporal coherence: gaussian falloff from centre
    tc = (
        (0.85 * np.exp(-r2 / 0.2) + 0.1 * RNG.random((NY, NX)))
        .clip(0, 1)
        .astype(np.float32)
    )

    recommended_mask = (tc >= 0.6).astype(np.float32)
    cc_labels = np.ones((NY, NX), dtype=np.float32)
    water_mask_arr = np.ones((NY, NX), dtype=np.float32)
    phase_sim = (0.9 * np.exp(-r2 / 0.15)).clip(0, 1).astype(np.float32)
    sw_disp = (displacement * 0.8).astype(np.float32)

    crs_wkt = _utm11n_crs_wkt()
    geotransform = f"{X0 - SPACING / 2} {SPACING} 0.0 {Y0 + SPACING / 2} 0.0 -{SPACING}"

    spatial_ref = xr.DataArray(
        0,
        attrs={
            "crs_wkt": crs_wkt,
            "semi_major_axis": 6378137.0,
            "semi_minor_axis": 6356752.314245179,
            "inverse_flattening": 298.257223563,
            "reference_ellipsoid_name": "WGS 84",
            "longitude_of_prime_meridian": 0.0,
            "prime_meridian_name": "Greenwich",
            "geographic_crs_name": "WGS 84",
            "horizontal_datum_name": "World Geodetic System 1984",
            "projected_crs_name": "WGS 84 / UTM zone 11N",
            "grid_mapping_name": "transverse_mercator",
            "latitude_of_projection_origin": 0.0,
            "longitude_of_central_meridian": -117.0,
            "false_easting": 500000.0,
            "false_northing": 0.0,
            "scale_factor_at_central_meridian": 0.9996,
            "GeoTransform": geotransform,
            "units": "unitless",
            "long_name": "Dummy variable with geo-referencing metadata in attributes",
        },
    )

    gm = {"grid_mapping": "spatial_ref"}
    ds = xr.Dataset(
        {
            "displacement": xr.DataArray(
                displacement,
                dims=["y", "x"],
                attrs={
                    "units": "meters",
                    "long_name": "Line-of-sight displacement",
                    "description": "Displacement along the radar LOS direction.",
                    **gm,
                },
            ),
            "short_wavelength_displacement": xr.DataArray(
                sw_disp,
                dims=["y", "x"],
                attrs={
                    "units": "meters",
                    "long_name": "Short wavelength displacement",
                    "wavelength_cutoff": "30000.0",
                    "wavelength_cutoff_units": "meters",
                    **gm,
                },
            ),
            "temporal_coherence": xr.DataArray(
                tc,
                dims=["y", "x"],
                attrs={"units": "unitless", "long_name": "Temporal Coherence", **gm},
            ),
            "recommended_mask": xr.DataArray(
                recommended_mask,
                dims=["y", "x"],
                attrs={
                    "units": "unitless",
                    "long_name": "Recommended Mask",
                    "temporal_coherence_threshold": "0.6",
                    **gm,
                },
            ),
            "connected_component_labels": xr.DataArray(
                cc_labels,
                dims=["y", "x"],
                attrs={
                    "units": "unitless",
                    "long_name": "Connected Component Labels",
                    **gm,
                },
            ),
            "water_mask": xr.DataArray(
                water_mask_arr,
                dims=["y", "x"],
                attrs={"units": "unitless", "long_name": "Water Mask", **gm},
            ),
            "phase_similarity": xr.DataArray(
                phase_sim,
                dims=["y", "x"],
                attrs={"units": "unitless", "long_name": "Phase Similarity", **gm},
            ),
            "estimated_phase_quality": xr.DataArray(
                phase_sim.copy(),
                dims=["y", "x"],
                attrs={
                    "units": "unitless",
                    "long_name": "Estimated phase quality",
                    **gm,
                },
            ),
            "persistent_scatterer_mask": xr.DataArray(
                np.zeros((NY, NX), dtype=np.float32),
                dims=["y", "x"],
                attrs={
                    "units": "unitless",
                    "long_name": "Persistent Scatterer Mask",
                    **gm,
                },
            ),
            "shp_counts": xr.DataArray(
                (20 * np.ones((NY, NX))).astype(np.float32),
                dims=["y", "x"],
                attrs={
                    "units": "unitless",
                    "long_name": "Statistically Homogeneous Pixels Counts",
                    **gm,
                },
            ),
            "timeseries_inversion_residuals": xr.DataArray(
                np.zeros((NY, NX), dtype=np.float32),
                dims=["y", "x"],
                attrs={
                    "units": "radians",
                    "long_name": "Timeseries inversion residuals",
                    **gm,
                },
            ),
            "displacement_corrected_constant_igs20": xr.DataArray(
                displacement.astype(np.float64),
                dims=["y", "x"],
                attrs={
                    "units": "meters",
                    "long_name": "Line-of-sight displacement",
                    **gm,
                },
            ),
            "reference_time": xr.DataArray(
                [np.datetime64(REF_DATETIME.replace(tzinfo=None), "ns")],
                dims=["time"],
                attrs={
                    "standard_name": "time",
                    "long_name": (
                        "Time corresponding to beginning of reference acquisition"
                    ),
                },
            ),
            "spatial_ref": spatial_ref,
        },
        coords={
            "x": xr.DataArray(x, dims=["x"], attrs={"units": "metres"}),
            "y": xr.DataArray(y, dims=["y"], attrs={"units": "metres"}),
            "time": [np.datetime64(SEC_DATETIME.replace(tzinfo=None), "ns")],
        },
        attrs={
            "Conventions": "CF-1.8",
            "contact": "opera-sds-ops@jpl.nasa.gov",
            "institution": "NASA JPL",
            "mission_name": "OPERA",
            "reference_document": "JPL D-108765",
            "title": "OPERA_L3_DISP-S1 Product",
        },
    )

    out = out_dir / DISP_FILENAME
    ds.to_netcdf(out, engine="h5netcdf")

    # Add /identification/radar_wavelength (Sentinel-1 C-band)
    import h5py

    with h5py.File(out, "a") as f:
        ident = f.create_group("identification")
        ident.create_dataset("radar_wavelength", data=np.float32(0.05546))
        ident["radar_wavelength"].attrs["units"] = "m"

    print(f"  created {out.name}")
    return out


# LOS GeoTIFF  (3 bands: east, north, up)
def create_los(out_dir: Path) -> Path:
    """Create a 3-band LOS ENU GeoTIFF and return its path."""
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import from_origin

    out_dir.mkdir(parents=True, exist_ok=True)

    # Descending Sentinel-1 geometry, typical for CA (incidence ~37°)
    los_up = np.full((NY, NX), 0.7986, dtype=np.float32)  # cos(37°)
    los_east = np.full((NY, NX), -0.5972, dtype=np.float32)  # ascending east
    los_north = np.full((NY, NX), 0.0744, dtype=np.float32)  # small north

    # Small spatial variation so the surface fitting has non-trivial geometry
    xx, yy = np.meshgrid(np.linspace(-0.02, 0.02, NX), np.linspace(-0.02, 0.02, NY))
    los_up += xx.astype(np.float32)
    los_east += yy.astype(np.float32)

    out = out_dir / LOS_FILENAME
    transform = from_origin(X0 - SPACING / 2, Y0 + SPACING / 2, SPACING, SPACING)
    crs = CRS.from_epsg(EPSG)

    with rasterio.open(
        out,
        "w",
        driver="GTiff",
        height=NY,
        width=NX,
        count=3,
        dtype=np.float32,
        crs=crs,
        transform=transform,
        nodata=0.0,
        compress="deflate",
    ) as dst:
        dst.write(los_east, 1)
        dst.write(los_north, 2)
        dst.write(los_up, 3)
        dst.update_tags(1, name="LOS East")
        dst.update_tags(2, name="LOS North")
        dst.update_tags(3, name="LOS Up")

    print(f"  created {out.name}")
    return out


# DEM GeoTIFF (static input)
def create_dem(out_dir: Path) -> Path:
    """Create a synthetic DEM GeoTIFF in UTM 11N and return its path.

    Uses the same grid as the synthetic DISP product (200 × 200, 30 m, UTM 11N)
    so the tropo pipeline needs no reprojection and shapes always match.
    """
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import from_origin

    out_dir.mkdir(parents=True, exist_ok=True)

    rng_dem = np.random.default_rng(7)
    xx, yy = np.meshgrid(np.linspace(0, 1, NX), np.linspace(0, 1, NY))
    elevation = (200 * xx + 150 * yy + 50 * rng_dem.random((NY, NX))).astype(np.float32)

    out = out_dir / DEM_FILENAME
    transform = from_origin(X0 - SPACING / 2, Y0 + SPACING / 2, SPACING, SPACING)
    crs = CRS.from_epsg(EPSG)

    with rasterio.open(
        out,
        "w",
        driver="GTiff",
        height=NY,
        width=NX,
        count=1,
        dtype=np.float32,
        crs=crs,
        transform=transform,
        nodata=-9999.0,
        compress="deflate",
    ) as dst:
        dst.write(elevation, 1)

    print(f"  created {out.name}")
    return out


# TROPO-ZENITH NetCDF files (reference + secondary)
def create_tropo(out_dir: Path) -> tuple[Path, Path]:
    """Create synthetic TROPO-ZENITH files for reference and secondary dates."""
    out_dir.mkdir(parents=True, exist_ok=True)

    n_height, n_lat, n_lon = 20, 50, 50
    height = np.linspace(0, 15_000, n_height)
    lat = np.linspace(34.5, 33.5, n_lat)  # north → south (decreasing)
    lon = np.linspace(-119.0, -117.0, n_lon)

    spatial_ref = xr.DataArray(
        0,
        attrs={
            "crs_wkt": (
                'GEOGCS["WGS 84",DATUM["WGS_1984",'
                'SPHEROID["WGS 84",6378137,298.257223563]],'
                'PRIMEM["Greenwich",0],'
                'UNIT["degree",0.0174532925199433]]'
            )
        },
    )

    rng_tropo = np.random.default_rng(99)
    paths = []

    for sensing_dt, filename in [
        (REF_DATETIME, REF_TROPO_FILENAME),
        (SEC_DATETIME, SEC_TROPO_FILENAME),
    ]:
        wet = (0.05 + 0.02 * rng_tropo.random((n_height, n_lat, n_lon))).astype(
            np.float32
        )
        hydro = (2.0 + 0.1 * rng_tropo.random((n_height, n_lat, n_lon))).astype(
            np.float32
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
                "time": [np.datetime64(sensing_dt.replace(tzinfo=None), "ns")],
            },
            attrs={"units": "meters", "model": "HRES"},
        )

        out = out_dir / filename
        ds.to_netcdf(out, engine="h5netcdf")
        print(f"  created {out.name}")
        paths.append(out)

    return paths[0], paths[1]


# GNSS lookup + .tenv8 files
def create_gnss(out_dir: Path) -> tuple[Path, Path]:
    """Return (lookup_file, tenv8_dir)."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # Place 4 stations inside the DISP UTM extent with a 20 % inset
    inset_x = NX * SPACING * 0.2
    inset_y = NY * SPACING * 0.2
    x_min = X0 + inset_x
    x_max = X0 + (NX - 1) * SPACING - inset_x
    y_max = Y0 - inset_y
    y_min = Y0 - (NY - 1) * SPACING + inset_y

    utm_pts = [
        (x_min, y_max),
        (x_max, y_max),
        (x_min, y_min),
        (x_max, y_min),
    ]

    # Convert to WGS84 lon/lat for the lookup file
    wgs84_pts = [_utm_to_wgs84(e, n) for e, n in utm_pts]

    # Write lookup table (format: {id:06d}  {lon}  {lat})
    lookup = out_dir / "grid_latlon_lookup.txt"
    with lookup.open("w") as f:
        for i, (lon, lat) in enumerate(wgs84_pts, 1):
            f.write(f"{i:06d}  {lon:.6f}  {lat:.6f}\n")

    # Create matching .tenv8 files
    # Format: year.frac  east(mm)  north(mm)  up(mm)  sig_e  sig_n  sig_u  flag
    # Span 2014.0 – 2018.0 at ~weekly cadence
    years = np.arange(2014.0, 2018.01, 0.0192)  # ~weekly (~52 obs/year)
    n = len(years)

    rng_gnss = np.random.default_rng(42)  # separate seed for GNSS noise

    # Linear velocity typical of California GNSS (IGS20 reference frame)
    v_east, v_north, v_up = 3.0, 0.0, -0.5  # mm/yr
    sig_e, sig_n, sig_u = 1.0, 1.0, 3.0

    for i in range(1, 5):
        path = out_dir / f"{i:06d}_IGS20.tenv8"
        t0 = 2016.0  # reference epoch
        east = v_east * (years - t0) + rng_gnss.normal(0, 0.5, n)
        north = v_north * (years - t0) + rng_gnss.normal(0, 0.5, n)
        up = v_up * (years - t0) + rng_gnss.normal(0, 1.5, n)

        with path.open("w") as f:
            for yr, e, no, u in zip(years, east, north, up):
                f.write(
                    f"{yr:.4f}  {e:11.3f}  {no:11.3f}  {u:11.3f}"
                    f"  {sig_e:.3f}  {sig_n:.3f}  {sig_u:.3f}  0\n"
                )

    print(f"  created {lookup.name} + 4 .tenv8 files")
    return lookup, out_dir


# Algorithm parameters YAML
def create_algorithm_params(out_dir: Path) -> Path:
    """Write the golden algorithm parameters YAML and return its path."""
    from cal_disp.config._algorithm import AlgorithmParameters, CalibrationOptions

    out_dir.mkdir(parents=True, exist_ok=True)
    params_file = out_dir / "algorithm_parameters.yaml"

    AlgorithmParameters(
        calibration_options=CalibrationOptions(
            posting_meters=SPACING,
            window_size_meters=15000.0,
            downsample_factor=1,
            calibration_surface_smoothing_sigma=0,
            grid_type="constant",
        )
    ).to_yaml(params_file, with_comments=False)

    print(f"  created {params_file}")
    return params_file


# Run workflow → golden output
def run_workflow(
    disp_file: Path,
    los_file: Path,
    dem_file: Path,
    gnss_lookup: Path,
    gnss_dir: Path,
    params_file: Path,
    ref_tropo_files: list[Path],
    sec_tropo_files: list[Path],
    output_dir: Path,
) -> Path:
    """Generate the golden output via the cal-disp CLI (config + run)."""
    import subprocess

    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "_work"
    work_dir.mkdir(parents=True, exist_ok=True)
    config_file = work_dir / "runconfig.yaml"

    # Build the tropo flags (one --ref/sec-tropo-files flag per file)
    ref_flags = [f for p in ref_tropo_files for f in ("--ref-tropo-files", str(p))]
    sec_flags = [f for p in sec_tropo_files for f in ("--sec-tropo-files", str(p))]

    subprocess.run(
        [
            "cal-disp",
            "config",
            "-d",
            str(disp_file),
            "-ul",
            str(gnss_lookup),
            "-ud",
            str(gnss_dir),
            "-uv",
            "0.3",
            "-ut",
            "constant",
            "--los-file",
            str(los_file),
            "--dem-file",
            str(dem_file),
            *ref_flags,
            *sec_flags,
            "-a",
            str(params_file),
            "--frame-id",
            str(FRAME_ID),
            "-o",
            str(output_dir),
            "--work-dir",
            str(work_dir),
            "-c",
            str(config_file),
        ],
        check=True,
    )

    subprocess.run(["cal-disp", "run", str(config_file)], check=True)

    output_files = sorted(output_dir.glob("OPERA_L4_DISP-CAL-S1_*.nc"))
    if not output_files:
        raise RuntimeError(f"No CalProduct found in {output_dir}")

    out = output_files[0]
    print(f"  golden output: {out}")
    return out


# Main
def main() -> None:
    """Generate all golden dataset files and run the calibration workflow."""
    global INPUT_DATA_DIR, CONFIGS_DIR, GOLDEN_OUTPUT_DIR

    args = _parse_args()
    root = _resolve_output_dir(args.output_dir) / "golden_datasets"

    INPUT_DATA_DIR = root / "input_data"
    CONFIGS_DIR = root / "configs"
    GOLDEN_OUTPUT_DIR = root / "golden_output"
    output_dir = root / "output"

    print("=== Generating golden dataset ===")
    print(f"  input_data    {INPUT_DATA_DIR}")
    print(f"  configs       {CONFIGS_DIR}")
    print(f"  golden_output {GOLDEN_OUTPUT_DIR}")
    print(f"  output        {output_dir}")
    print()

    print("Creating DISP product...")
    disp = create_disp(INPUT_DATA_DIR / "disp")

    print("Creating LOS GeoTIFF...")
    los = create_los(INPUT_DATA_DIR / "static_input")

    print("Creating DEM GeoTIFF...")
    dem = create_dem(INPUT_DATA_DIR / "static_input")

    print("Creating TROPO files...")
    ref_tropo, sec_tropo = create_tropo(INPUT_DATA_DIR / "tropo")

    print("Creating GNSS files...")
    gnss_lookup, gnss_dir = create_gnss(INPUT_DATA_DIR / "gnss")

    print("Creating algorithm parameters...")
    params = create_algorithm_params(CONFIGS_DIR)

    print("Creating output directory...")
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Running calibration workflow...")
    try:
        run_workflow(
            disp,
            los,
            dem,
            gnss_lookup,
            gnss_dir,
            params,
            [ref_tropo],
            [sec_tropo],
            GOLDEN_OUTPUT_DIR,
        )
    except Exception as exc:
        print(f"\nERROR running workflow: {exc}", file=sys.stderr)
        print(
            "Input files were created successfully.  Fix the error above and"
            " re-run to generate the golden output.",
            file=sys.stderr,
        )
        sys.exit(1)

    print()
    print(f"Done.  Golden data generated at {root}")


if __name__ == "__main__":
    main()
