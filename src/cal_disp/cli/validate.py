from __future__ import annotations

from pathlib import Path

import click


@click.command("validate")
@click.argument("reference", type=click.Path(exists=True, path_type=Path))
@click.argument("test", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--tolerance",
    "-t",
    type=float,
    default=1e-6,
    show_default=True,
    help="Tolerance for floating point comparison.",
)
@click.option(
    "--group",
    "-g",
    type=click.Choice(["main", "auxiliary", "all"]),
    default="all",
    show_default=True,
    help="Which group to validate: main (required data), auxiliary (optional), or all.",
)
@click.pass_context
def validate_cli(
    ctx: click.Context,
    reference: Path,
    test: Path,
    tolerance: float,
    group: str,
) -> None:
    """Validate a DISP-CAL product against a reference.

    Compares TEST product against REFERENCE product and reports differences.
    Validates structure (dimensions, types) and data values within tolerance.

    Examples
    --------
    Basic validation:

        cal-disp validate reference.nc test.nc

    With custom tolerance:

        cal-disp validate reference.nc test.nc --tolerance 1e-5

    Validate only main group:

        cal-disp validate reference.nc test.nc --group main

    """
    from cal_disp._log import setup_logging
    from cal_disp.validate import compare_cal_products

    debug = ctx.obj.get("debug", False)
    setup_logging(logger_name="cal_disp", level="DEBUG" if debug else "INFO")

    success = compare_cal_products(
        reference_file=reference,
        test_file=test,
        tolerance=tolerance,
        group=group,
    )

    if not success:
        raise click.Abort()


@click.command("validate-golden")
@click.argument(
    "golden_dir", type=click.Path(exists=True, file_okay=False, path_type=Path)
)
@click.option(
    "--tolerance",
    "-t",
    type=float,
    default=1e-6,
    show_default=True,
    help="Tolerance for floating point comparison.",
)
@click.option(
    "--group",
    "-g",
    type=click.Choice(["main", "auxiliary", "all"]),
    default="all",
    show_default=True,
    help="Which group to validate: main, auxiliary, or all.",
)
@click.pass_context
def validate_golden_cli(
    ctx: click.Context,
    golden_dir: Path,
    tolerance: float,
    group: str,
) -> None:
    r"""Re-run the calibration on golden inputs and compare against the reference.

    GOLDEN_DIR is the root directory produced by build_golden_output.sh.
    It must contain:

    \b
      input_data/    calibration inputs (disp, gnss, static_input, tropo)
      configs/       algorithm_parameters.yaml
      golden_output/ known-good reference CalProduct .nc
      output/        populated during this run

    Examples
    --------
    Validate a delivered golden dataset:

        cal-disp validate-golden /path/to/golden_dataset

    With custom tolerance:

        cal-disp validate-golden /path/to/golden_dataset --tolerance 1e-5

    """
    import subprocess

    from cal_disp._log import setup_logging

    debug = ctx.obj.get("debug", False)
    setup_logging(logger_name="cal_disp", level="DEBUG" if debug else "INFO")

    golden_dir = golden_dir.resolve()
    input_dir = golden_dir / "input_data"
    configs_dir = golden_dir / "configs"
    reference_dir = golden_dir / "golden_output"
    output_dir = golden_dir / "output"
    work_dir = output_dir / "_work"

    # Validate directory layout
    for d in (input_dir, configs_dir, reference_dir):
        if not d.is_dir():
            raise click.ClickException(
                f"Expected directory not found: {d}\n"
                "Make sure GOLDEN_DIR was produced by build_golden_output.sh."
            )

    import re

    # Locate required input files
    disp_files = sorted((input_dir / "disp").glob("OPERA_L3_DISP-S1_*.nc"))
    los_files = sorted((input_dir / "static_input").glob("*line_of_sight_enu.tif"))
    dem_files = sorted((input_dir / "static_input").glob("*_dem.tif"))
    tropo_files = sorted((input_dir / "tropo").glob("OPERA_L4_TROPO-ZENITH_*.nc"))
    algo_file = configs_dir / "algorithm_parameters.yaml"
    ref_files = sorted(reference_dir.glob("OPERA_L4_DISP-CAL-S1_*.nc"))

    # Lookup file encodes its version: grid_latlon_lookup_v{ver}.txt
    gnss_dir = input_dir / "gnss"
    lookup_files = sorted(gnss_dir.glob("grid_latlon_lookup_v*.txt"))
    lookup_file = lookup_files[0] if lookup_files else Path()

    for label, collection in [
        ("DISP file", disp_files),
        ("LOS file", los_files),
        ("DEM file", dem_files),
        ("TROPO files (need 2)", tropo_files if len(tropo_files) == 2 else []),
        (
            "UNR lookup (grid_latlon_lookup_v*.txt)",
            [lookup_file] if lookup_file.exists() else [],
        ),
        ("algorithm params", [algo_file] if algo_file.exists() else []),
        ("golden reference", ref_files),
    ]:
        if not collection:
            raise click.ClickException(f"Missing {label} in {golden_dir}")

    disp_file = disp_files[0]
    frame_id = disp_file.stem.split("_F")[1].split("_")[0].lstrip("0") or "0"

    # Parse UNR version from lookup filename (e.g. grid_latlon_lookup_v0.3.txt → "0.3")
    _v_match = re.search(r"v(\d+(?:\.\d+)+)", lookup_file.name)
    if not _v_match:
        raise click.ClickException(
            f"Cannot determine UNR version from lookup filename '{lookup_file.name}'. "
            "Expected format: grid_latlon_lookup_v<version>.txt"
        )
    unr_version = _v_match.group(1)

    # Parse UNR type from algorithm parameters
    unr_type = "constant"
    algo_text = algo_file.read_text()
    for line in algo_text.splitlines():
        if "grid_type" in line:
            unr_type = line.split(":")[-1].strip().strip('"').strip("'")
            break

    work_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    config_file = work_dir / "runconfig.yaml"

    ref_flags = [f for p in [tropo_files[0]] for f in ("--ref-tropo-files", str(p))]
    sec_flags = [f for p in [tropo_files[1]] for f in ("--sec-tropo-files", str(p))]

    click.echo("[1/3] Generating configuration...")
    subprocess.run(
        [
            "cal-disp",
            "config",
            "-d",
            str(disp_file),
            "-ul",
            str(lookup_file),
            "-ud",
            str(input_dir / "gnss"),
            "-uv",
            unr_version,
            "-ut",
            unr_type,
            "--los-file",
            str(los_files[0]),
            "--dem-file",
            str(dem_files[0]),
            *ref_flags,
            *sec_flags,
            "-a",
            str(algo_file),
            "--frame-id",
            frame_id,
            "-o",
            str(output_dir),
            "--work-dir",
            str(work_dir),
            "-c",
            str(config_file),
        ],
        check=True,
    )

    click.echo("[2/3] Running calibration...")
    subprocess.run(
        ["cal-disp", "run", str(config_file)],
        check=True,
    )

    test_files = sorted(output_dir.glob("OPERA_L4_DISP-CAL-S1_*.nc"))
    if not test_files:
        raise click.ClickException("Calibration produced no output NetCDF.")

    click.echo("[3/3] Validating output...")
    from cal_disp.validate import compare_cal_products

    success = compare_cal_products(
        reference_file=ref_files[0],
        test_file=test_files[0],
        tolerance=tolerance,
        group=group,
    )

    # Clean up scratch directory
    import shutil

    shutil.rmtree(work_dir, ignore_errors=True)

    if success:
        click.echo("\nValidation passed.")
    else:
        click.echo("\nValidation FAILED.", err=True)
        raise click.Abort()
