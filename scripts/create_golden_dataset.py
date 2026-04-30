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

Outputs (written under <output-dir>/)
--------------------------------------
golden/
    disp/          -- 200 × 200 synthetic DISP-S1 NetCDF (30 m, UTM 11N)
    gnss/          -- UNR lookup table + 4 × .tenv8 files
    los.tif        -- 3-band LOS GeoTIFF
    water_mask.tif
    algorithm_parameters.yaml

golden_output/
    OPERA_L4_CAL-DISP-S1_IW_F36540_VV_*.nc  -- expected CalProduct

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
GOLDEN_DATA_DIR: Path
GOLDEN_OUTPUT_DIR: Path

# Grid constants — match real frame-36540 resolution and CRS
NY, NX = 200, 200
SPACING = 30.0  # metres (same as real data)
EPSG = 32611  # WGS 84 / UTM zone 11N

# Top-left pixel *centre* coords (UTM 11N metres)
# Chosen so that synthetic GNSS stations sit inside the extent.
X0 = 405000.0
Y0 = 3778000.0

# Acquisition dates taken directly from real frame-36540 test file
REF_DATETIME = datetime(2016, 7, 24, 1, 58, 9, tzinfo=timezone.utc)
SEC_DATETIME = datetime(2016, 8, 5, 1, 58, 9, tzinfo=timezone.utc)
FRAME_ID = 36540
DISP_FILENAME = (
    f"OPERA_L3_DISP-S1_IW_F{FRAME_ID:05d}_VV_"
    f"{REF_DATETIME:%Y%m%dT%H%M%S}Z_"
    f"{SEC_DATETIME:%Y%m%dT%H%M%S}Z_"
    "v1.0_20250101T000000Z.nc"
)

RNG = np.random.default_rng(0)  # fully deterministic


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
    ramp = 0.005 * xx + 0.003 * yy  # metres

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

    out = out_dir / "los.tif"
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


# Water mask GeoTIFF


def create_water_mask(out_dir: Path) -> Path:
    """Create an all-land water mask GeoTIFF and return its path."""
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import from_origin

    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "water_mask.tif"
    transform = from_origin(X0 - SPACING / 2, Y0 + SPACING / 2, SPACING, SPACING)
    crs = CRS.from_epsg(EPSG)
    mask = np.ones((NY, NX), dtype=np.uint8)  # all land

    with rasterio.open(
        out,
        "w",
        driver="GTiff",
        height=NY,
        width=NX,
        count=1,
        dtype=np.uint8,
        crs=crs,
        transform=transform,
        nodata=255,
        compress="deflate",
    ) as dst:
        dst.write(mask, 1)

    print(f"  created {out.name}")
    return out


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
            window_size_meters=1500.0,  # 50-pixel window on 200×200 grid
            downsample_factor=1,
            calibration_surface_smoothing_sigma=0,
            grid_type="constant",  # deterministic; no epoch-specific GNSS call
        )
    ).to_yaml(params_file, with_comments=False)

    print(f"  created {params_file}")
    return params_file


# Run workflow → golden output


def run_workflow(
    disp_file: Path,
    los_file: Path,
    gnss_lookup: Path,
    gnss_dir: Path,
    params_file: Path,
    output_dir: Path,
) -> Path:
    """Run the calibration workflow and return the path to the output CalProduct."""
    from cal_disp.config._algorithm import AlgorithmParameters
    from cal_disp.workflow import run_calibration

    output_dir.mkdir(parents=True, exist_ok=True)

    params = AlgorithmParameters.from_yaml(params_file)

    out = run_calibration(
        disp_file=disp_file,
        unr_grid_latlon_file=gnss_lookup,
        unr_timeseries_dir=gnss_dir,
        output_dir=output_dir,
        algorithm_parameters=params,
        los_file=los_file,
        work_directory=output_dir / "_work",
        calibration_reference_type="constant",
        calibration_reference_reference_frame="IGS20",
    )
    print(f"  golden output: {out}")
    return out


# Main


def main() -> None:
    """Generate all golden dataset files and run the calibration workflow."""
    global GOLDEN_DATA_DIR, GOLDEN_OUTPUT_DIR

    args = _parse_args()
    root = _resolve_output_dir(args.output_dir)

    GOLDEN_DATA_DIR = root / "golden"
    GOLDEN_OUTPUT_DIR = root / "golden_output"

    print("=== Generating golden dataset ===")
    print(f"  inputs  → {GOLDEN_DATA_DIR}")
    print(f"  output  → {GOLDEN_OUTPUT_DIR}")
    print()

    print("Creating DISP product...")
    disp = create_disp(GOLDEN_DATA_DIR / "disp")

    print("Creating LOS GeoTIFF...")
    los = create_los(GOLDEN_DATA_DIR)

    print("Creating water mask...")
    create_water_mask(GOLDEN_DATA_DIR)

    print("Creating GNSS files...")
    gnss_lookup, gnss_dir = create_gnss(GOLDEN_DATA_DIR / "gnss")

    print("Creating algorithm parameters...")
    params = create_algorithm_params(GOLDEN_DATA_DIR)

    print("Running calibration workflow...")
    try:
        run_workflow(disp, los, gnss_lookup, gnss_dir, params, GOLDEN_OUTPUT_DIR)
    except Exception as exc:
        print(f"\nERROR running workflow: {exc}", file=sys.stderr)
        print(
            "Input files were created successfully.  Fix the error above and"
            " re-run to generate the golden output.",
            file=sys.stderr,
        )
        sys.exit(1)

    print()
    print(f"Done.  Set ${_ENV_VAR}={root} to enable integration tests.")


if __name__ == "__main__":
    main()
