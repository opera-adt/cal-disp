from __future__ import annotations

from pathlib import Path

from cal_disp.product import (
    DispProduct,
    StaticLayer,
    TropoProduct,
    compute_los_correction,
    interpolate_in_time,
    interpolate_to_dem_surface,
)


def prepare_troposphere_correction(
    disp_file: Path,
    dem_file: Path,
    los_file: Path,
    reference_tropo_files: list[Path],
    secondary_tropo_files: list[Path],
    output_dir: Path,
    bounds_buffer: float = 0.2,
    max_time_offset_hours: float = 6.0,
) -> list[Path]:
    """Prepare tropospheric delay corrections for displacement products.

    Loads DEM and line-of-sight (LOS) data, interpolates tropospheric delay
    to the target dates, projects to LOS, and saves the results.

    Parameters
    ----------
    disp_file : Path
         Path to OPERA Displacement product.
    dem_file : Path
        Path to digital elevation model file.
    los_file : Path
        Path to line-of-sight geometry file.
    reference_tropo_files : list[Path]
        Troposphere files for primary date (1-2 files for temporal interpolation).
    secondary_tropo_files : list[Path]
        Troposphere files for secondary date (1-2 files for temporal interpolation).
    output_dir : Path
        Directory where corrected troposphere layers will be saved.
    bounds_buffer : float, optional
        Buffer factor for spatial bounds query. Default is 0.2.
    max_time_offset_hours : float, optional
        Maximum allowed time offset between tropo file and target date. Default is 6.0.

    Returns
    -------
    tuple[Path, Path]
        Paths to the saved reference and secondary troposphere corrections.

    Raises
    ------
    ValueError
        If troposphere files don't match target dates within time threshold,
        or if more than 2 files are provided per date.

    Notes
    -----
    Troposphere files are temporally interpolated if two files are provided
    (one before and one after the target datetime). If only one file is provided,
    it is used directly without temporal interpolation.

    Examples
    --------
    >>> ref_files = [Path("tropo_20240101_00.nc"), Path("tropo_20240101_12.nc")]
    >>> sec_files = [Path("tropo_20240107_00.nc")]
    >>> ref_path, sec_path = prepare_troposphere_correction(
    ...     disp_product=my_product,
    ...     dem_file=Path("dem.tif"),
    ...     los_file=Path("los.tif"),
    ...     reference_tropo_files=ref_files,
    ...     secondary_tropo_files=sec_files,
    ...     output_dir=Path("corrections/troposphere"),
    ... )

    """
    # LOAD DISP product
    disp_product = DispProduct.from_path(disp_file)

    # Load DEM and add CRS
    dem = StaticLayer.from_path(dem_file).to_dataset()
    dem = dem.rio.write_crs(dem.attrs["crs_wkt"])
    max_height = dem.dem.max().values + 3e3

    # Load LOS geometry
    los_ds = StaticLayer.from_path(los_file).to_dataset()

    # Get spatial bounds from product
    bounds = list(disp_product.get_bounds_wgs84().values())

    # Create output directory
    tropo_dir = output_dir / "troposphere_delay"
    tropo_dir.mkdir(exist_ok=True, parents=True)

    output_paths = []

    # Process both reference and secondary dates
    for tropo_files, target_date in [
        (reference_tropo_files, disp_product.primary_date),
        (secondary_tropo_files, disp_product.secondary_date),
    ]:
        # Validate temporal matching
        for tropo_file in tropo_files:
            tropo_prod = TropoProduct.from_path(tropo_file)
            if not tropo_prod.matches_date(target_date, hours=max_time_offset_hours):
                raise ValueError(
                    f"{tropo_file.stem} is not within {max_time_offset_hours}h "
                    f"of target date {target_date}"
                )

        # Get troposphere delay (with temporal interpolation if needed)
        if len(tropo_files) == 2:
            tropo_early, tropo_late = tropo_files
            tropo_ds = interpolate_in_time(
                TropoProduct.from_path(tropo_early),
                TropoProduct.from_path(tropo_late),
                target_datetime=target_date,
                bounds=bounds,
                max_height=max_height,
                bounds_buffer=bounds_buffer,
            )
        elif len(tropo_files) == 1:
            tropo_ds = TropoProduct.from_path(tropo_files[0]).get_total_delay(
                bounds=bounds,
                max_height=max_height,
                bounds_buffer=bounds_buffer,
            )
        else:
            raise ValueError(
                f"Expected 1-2 troposphere files per date, got {len(tropo_files)}"
            )

        # Interpolate to DEM surface
        tropo_dem = interpolate_to_dem_surface(tropo_ds, dem["dem"])

        # Project to LOS and save
        output_path = tropo_dir / f"{target_date.strftime('%Y%m%d')}.tif"
        compute_los_correction(
            tropo_dem,
            los_ds["los_up"],
            output_path=output_path,
        )
        output_paths.append(output_path)

    return output_paths
