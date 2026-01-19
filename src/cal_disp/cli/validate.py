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
    """Validate a CAL-DISP product against a reference.

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
