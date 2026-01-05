from __future__ import annotations

import click

from .config import config_cli
from .download import download_group
from .run import run_cli
from .validate import validate_cli


@click.group(name="cal-disp")
@click.version_option()
@click.option("--debug", is_flag=True, help="Enable debug logging.")
@click.pass_context
def cli_app(ctx: click.Context, debug: bool) -> None:
    """OPERA-DISP calibration workflow toolkit."""
    ctx.ensure_object(dict)
    ctx.obj["debug"] = debug


cli_app.add_command(config_cli)
cli_app.add_command(run_cli)
cli_app.add_command(validate_cli)
cli_app.add_command(download_group)


if __name__ == "__main__":
    cli_app()
