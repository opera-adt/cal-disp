from __future__ import annotations

from pathlib import Path

import click

__all__ = ["run_cli", "run_main"]


def run_main(config_file: str | Path, debug: bool = False) -> Path:
    """Run the calibration workflow.

    Parameters
    ----------
    config_file : str or Path
        Path to YAML configuration file.
    debug : bool, optional
        Enable debug logging. Default is False.

    Returns
    -------
    Path
        Path to output calibrated displacement file.

    """
    # Lazy imports so --help is fast
    from cal_disp.config.pge_runconfig import RunConfig
    from cal_disp.main import run

    pge_runconfig = RunConfig.from_yaml_file(Path(config_file))
    runconfig = pge_runconfig.to_workflow()
    return run(runconfig=runconfig, debug=debug)


@click.command("run")
@click.argument("config_file", type=click.Path(exists=True, path_type=Path))
@click.pass_context
def run_cli(ctx: click.Context, config_file: Path) -> None:
    """Run the calibration workflow for CONFIG_FILE.

    CONFIG_FILE is a YAML configuration file specifying input files,
    output paths, and processing parameters.

    """
    run_main(config_file=config_file, debug=ctx.obj["debug"])
