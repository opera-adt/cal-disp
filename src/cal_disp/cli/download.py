from __future__ import annotations

from datetime import datetime
from pathlib import Path

import click


@click.group(name="download")
def download_group():
    """Sub-commands for downloading prerequisite data."""
    from cal_disp._log import setup_logging

    setup_logging(logger_name="cal_disp")


@download_group.command()
@click.option(
    "--frame-id",
    type=int,
    required=True,
    help="OPERA DISP-S1 frame ID.",
    show_default=True,
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Directory to save downloaded DISP products.",
)
@click.option(
    "--start",
    "-s",
    type=click.DateTime(formats=["%Y-%m-%d"]),
    default=None,
    help="Start date (YYYY-MM-DD). Based on secondary date of DISP.",
)
@click.option(
    "--end",
    "-e",
    type=click.DateTime(formats=["%Y-%m-%d"]),
    default=None,
    help="End date (YYYY-MM-DD). Based on secondary date of DISP.",
)
@click.option(
    "--num-workers",
    "-n",
    type=int,
    default=1,
    help="Number of parallel download workers.",
    show_default=True,
)
def disp_s1(
    frame_id: int,
    output_dir: Path,
    start: datetime | None,
    end: datetime | None,
    num_workers: int,
) -> None:
    r"""Download OPERA DISP-S1 products for a frame.

    Downloads displacement products from the OPERA DISP-S1 archive for
    the specified frame and date range. Products are filtered based on
    the secondary date of each interferogram.

    Examples
    --------
    Download all products for a frame:
        $ cal-disp download disp-s1 --frame-id 8882 -o ./disp_data

    Download products for specific date range:
        $ cal-disp download disp-s1 --frame-id 8882 -o ./disp_data \\
            --start 2022-07-01 --end 2022-07-31

    Use more workers for faster downloads:
        $ cal-disp download disp-s1 --frame-id 8882 -o ./disp_data -n 8

    """
    from cal_disp.download import download_disp

    download_disp(
        frame_id=frame_id,
        output_dir=output_dir,
        start=start,
        end=end,
        num_workers=num_workers,
    )
    click.echo(f"Download complete: {output_dir}")


@download_group.command()
@click.option(
    "--frame-id",
    type=int,
    required=True,
    help="OPERA DISP-S1 frame ID.",
    show_default=True,
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Directory to save downloaded UNR parquet file.",
)
@click.option(
    "--margin",
    "-m",
    type=float,
    default=0.5,
    help="Margin in degrees to expand frame bounding box.",
    show_default=True,
)
def unr(
    frame_id: int,
    output_dir: Path,
    margin: float,
) -> None:
    r"""Download UNR GPS timeseries data for a DISP-S1 frame.

    Downloads GPS timeseries grid from the Nevada Geodetic Laboratory (UNR)
    within the frame's bounding box (expanded by margin). Data is saved
    as a parquet file for efficient loading. Output directory is created
    if it doesn't exist.

    Examples
    --------
    Download UNR data for a frame:
        $ cal-disp download unr --frame-id 8882 -o ./unr_data

    Expand bounding box by 1 degree:
        $ cal-disp download unr --frame-id 8882 -o ./unr_data -m 1.0

    """
    from cal_disp.download import download_unr_grid

    output_dir.mkdir(exist_ok=True, parents=True)

    download_unr_grid(
        frame_id=frame_id,
        output_dir=output_dir,
        margin_deg=margin,
    )
    click.echo(f"Download complete: {output_dir}")


@download_group.command()
@click.option(
    "--input-file",
    "-i",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Path to OPERA DISP-S1 product (.nc file).",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Directory to save downloaded TROPO products.",
)
@click.option(
    "--num-workers",
    "-n",
    type=int,
    default=4,
    help="Number of parallel download workers.",
    show_default=True,
)
@click.option(
    "--interp",
    is_flag=True,
    help="Download 2 scenes per date for temporal interpolation.",
)
def tropo(
    input_file: Path,
    output_dir: Path,
    num_workers: int,
    interp: bool,
) -> None:
    """Download OPERA TROPO for a DISP-S1 file.

    Extracts sensing times from the input DISP-S1 product filename and
    downloads corresponding OPERA TROPO data. Output directory
    is created if it doesn't exist.

    The input file must follow OPERA naming convention:
    OPERA_L3_DISP-S1_IW_F{frame}_VV_{ref_date}_{sec_date}_v1.0_{proc_date}.nc

    Examples
    --------
    Basic usage:
        $ cal-disp download tropo -i disp_product.nc -o ./tropo_data

    With temporal interpolation:
        $ cal-disp download tropo -i disp_product.nc -o ./tropo_data --interp

    Using more workers:
        $ cal-disp download tropo -i disp_product.nc -o ./tropo_data -n 8

    """
    from cal_disp.download import download_tropo
    from cal_disp.download.utils import extract_sensing_times_from_file

    sensing_times = extract_sensing_times_from_file(input_file)
    output_dir.mkdir(exist_ok=True, parents=True)

    download_tropo(
        disp_times=sensing_times,
        output_dir=output_dir,
        num_workers=num_workers,
        interp=interp,
    )
    click.echo(f"Download complete: {output_dir}")


@download_group.command(name="burst-bounds")
@click.option(
    "--input-file",
    "-i",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    required=True,
    help="Path to OPERA DISP-S1 product (.nc file).",
)
@click.option(
    "--output-dir",
    "-o",
    type=click.Path(file_okay=False, path_type=Path),
    required=True,
    help="Directory to save burst boundary tiles.",
)
def burst_bounds(
    input_file: Path,
    output_dir: Path,
) -> None:
    r"""Download S1 CSLC boundary tiles for a DISP-S1 file.

    Extracts frame ID and sensing times from the DISP-S1 filename,
    then downloads corresponding Sentinel-1 burst boundary geometries.
    Output directory is created if it doesn't exist.

    Only supports DISP-S1 products (sensor must be S1).

    Examples
    --------
    Basic usage:
        $ cal-disp download burst-bounds -i disp_product.nc -o ./burst_data

    Full path example:
        $ cal-disp download burst-bounds \\
        -i OPERA_L3_DISP-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0.nc \\
        -o ./burst_bounds

    """
    from cal_disp.download import generate_s1_burst_tiles
    from cal_disp.download.utils import extract_sensing_times_from_file

    # Parse filename: OPERA_L3_DISP-S1_IW_F{frame}_VV_{dates}...
    parts = input_file.stem.split("_")
    sensor = parts[2].split("-")[1]  # DISP-S1 -> S1
    frame_id = int(parts[4].lstrip("F"))  # F08882 -> 8882

    if sensor != "S1":
        raise click.ClickException(f"Only DISP-S1 products supported, got: {sensor}")

    sensing_times = extract_sensing_times_from_file(input_file)
    output_dir.mkdir(exist_ok=True, parents=True)

    for sensing_time in sensing_times:
        click.echo(f"Downloading burst bounds for {sensing_time}")
        generate_s1_burst_tiles(
            frame_id=frame_id,
            sensing_time=sensing_time,
            output_dir=output_dir,
        )

    click.echo(f"Download complete: {output_dir}")
