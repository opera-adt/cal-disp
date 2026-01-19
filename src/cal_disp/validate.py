from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import xarray as xr

from cal_disp.product import CalProduct

logger = logging.getLogger(__name__)


def compare_cal_products(
    reference_file: Path,
    test_file: Path,
    tolerance: float = 1e-6,
    group: str = "all",
) -> bool:
    """Compare two CAL-DISP products.

    Parameters
    ----------
    reference_file : Path
        Reference (golden) product file.
    test_file : Path
        Test product file to validate.
    tolerance : float, optional
        Tolerance for floating point comparison. Default is 1e-6.
    group : str, optional
        Which group to validate: "main", "auxiliary", or "all". Default is "all".

    Returns
    -------
    bool
        True if products match within tolerance, False otherwise.

    """
    logger.info(f"Comparing {test_file.name} against {reference_file.name}")

    # Load products
    ref_cal = CalProduct.from_path(reference_file)
    test_cal = CalProduct.from_path(test_file)

    # Compare metadata
    if not _compare_metadata(ref_cal, test_cal):
        return False

    # Validate structure (dimensions, list conversions)
    if not _validate_product_structure(test_cal):
        return False

    # Compare main group
    if group in ("main", "all"):
        if not _compare_group(ref_cal, test_cal, "main", tolerance):
            return False

    # Compare model_3d group if it exists
    if group in ("auxiliary", "all"):
        ref_has_model = ref_cal.has_auxiliary()
        test_has_model = test_cal.has_auxiliary()

        if ref_has_model != test_has_model:
            logger.error(
                f"auxiliary group mismatch: reference={ref_has_model},"
                f" test={test_has_model}"
            )
            return False

        if ref_has_model and not _compare_group(
            ref_cal, test_cal, "auxiliary", tolerance
        ):
            return False

    logger.info("✓ Validation passed")
    return True


def _validate_product_structure(cal: CalProduct) -> bool:
    """Validate product structure (dimensions, types).

    Ensures:
    - Calibration has (time, y, x) dimensions
    - List variables are strings, not arrays
    - No unnamed dimensions (dim_0)
    """
    logger.info("Validating product structure...")

    ds = cal.open_dataset()
    all_valid = True

    # Check calibration dimensions
    if "calibration" in ds.data_vars:
        cal_dims = ds["calibration"].dims
        if cal_dims != ("time", "y", "x"):
            logger.error(
                f"calibration has wrong dimensions: {cal_dims}, expected (time, y, x)"
            )
            all_valid = False
        else:
            logger.info("  ✓ calibration dimensions")

    # Check for unnamed dimensions
    for var in ds.data_vars:
        if any("dim_" in str(dim) for dim in ds[var].dims):
            logger.error(f"Variable '{var}' has unnamed dimension: {ds[var].dims}")
            all_valid = False

    # Check identification group if it exists
    try:
        with xr.open_dataset(cal.path, group="identification") as id_ds:
            if not _validate_identification_structure(id_ds):
                all_valid = False
    except (KeyError, OSError):
        # Identification group might not exist in all products
        pass

    return all_valid


def _validate_identification_structure(ds: xr.Dataset) -> bool:
    """Validate identification group structure.

    Ensures list variables are scalar strings, not arrays.
    """
    list_vars = [
        "source_calibration_file_list",
        "source_data_file_list",
        "source_data_satellite_names",
    ]

    all_valid = True
    for var in list_vars:
        if var not in ds.data_vars:
            continue

        data_arr = ds[var]

        # Should be scalar (no dimensions)
        if data_arr.dims != ():
            logger.error(
                f"'{var}' should be scalar but has dimensions: {data_arr.dims}"
            )
            all_valid = False
            continue

        # Should be string type
        value = data_arr.item()
        if not isinstance(value, str):
            logger.error(f"'{var}' should be string but is {type(value)}")
            all_valid = False
            continue

        logger.info(f"  ✓ {var} is scalar string")

    return all_valid


def _compare_metadata(ref: CalProduct, test: CalProduct) -> bool:
    """Compare product metadata."""
    checks = {
        "frame_id": (ref.frame_id, test.frame_id),
        "sensor": (ref.sensor, test.sensor),
        "mode": (ref.mode, test.mode),
        "polarization": (ref.polarization, test.polarization),
        "reference_date": (ref.reference_date, test.reference_date),
        "secondary_date": (ref.secondary_date, test.secondary_date),
        "version": (ref.version, test.version),
    }

    all_match = True
    for field, (ref_val, test_val) in checks.items():
        if ref_val != test_val:
            logger.error(
                f"Metadata mismatch in {field}: reference={ref_val}, test={test_val}"
            )
            all_match = False

    return all_match


def _compare_group(
    ref: CalProduct,
    test: CalProduct,
    group: str,
    tolerance: float,
) -> bool:
    """Compare a specific group."""
    group_name = "main" if group == "main" else "auxiliary"
    logger.info(f"Validating {group_name} group...")

    # Open datasets
    if group == "main":
        ref_ds = ref.open_dataset()
        test_ds = test.open_dataset()
    else:
        ref_ds = ref.open_auxiliary()
        test_ds = test.open_auxiliary()

    # Compare coordinates
    if not _compare_coords(ref_ds, test_ds):
        return False

    # Compare data variables
    ref_vars = set(ref_ds.data_vars) - {"spatial_ref"}
    test_vars = set(test_ds.data_vars) - {"spatial_ref"}

    if ref_vars != test_vars:
        logger.error(f"Variable mismatch in {group_name}:")
        logger.error(f"  Reference only: {ref_vars - test_vars}")
        logger.error(f"  Test only: {test_vars - ref_vars}")
        return False

    # Compare each variable
    all_match = True
    for var in ref_vars:
        if not _compare_array(ref_ds[var], test_ds[var], var, tolerance):
            all_match = False

    return all_match


def _compare_coords(ref_ds: xr.Dataset, test_ds: xr.Dataset) -> bool:
    """Compare dataset coordinates."""
    if set(ref_ds.coords) != set(test_ds.coords):
        logger.error("Coordinate mismatch:")
        logger.error(f"  Reference: {set(ref_ds.coords)}")
        logger.error(f"  Test: {set(test_ds.coords)}")
        return False

    all_match = True
    for coord in ref_ds.coords:
        if coord == "spatial_ref":
            continue

        ref_vals = ref_ds[coord].values
        test_vals = test_ds[coord].values

        # Handle datetime coordinates separately
        if np.issubdtype(ref_vals.dtype, np.datetime64):
            if not np.array_equal(ref_vals, test_vals):
                logger.error(f"Coordinate '{coord}' values differ")
                all_match = False
        else:
            # Numeric coordinates
            if not np.allclose(ref_vals, test_vals, rtol=1e-9):
                logger.error(f"Coordinate '{coord}' values differ")
                all_match = False

    return all_match


def _compare_array(
    ref_arr: xr.DataArray,
    test_arr: xr.DataArray,
    name: str,
    tolerance: float,
) -> bool:
    """Compare two data arrays."""
    # Check dimensions match
    if ref_arr.dims != test_arr.dims:
        logger.error(
            f"Dimension mismatch for '{name}': {ref_arr.dims} vs {test_arr.dims}"
        )
        return False

    ref_data = ref_arr.values
    test_data = test_arr.values

    # Check shape
    if ref_data.shape != test_data.shape:
        logger.error(
            f"Shape mismatch for '{name}': {ref_data.shape} vs {test_data.shape}"
        )
        return False

    # Handle NaNs
    ref_nan_mask = np.isnan(ref_data)
    test_nan_mask = np.isnan(test_data)

    if not np.array_equal(ref_nan_mask, test_nan_mask):
        logger.error(f"NaN pattern differs for '{name}'")
        logger.error(f"  Reference NaNs: {ref_nan_mask.sum()}")
        logger.error(f"  Test NaNs: {test_nan_mask.sum()}")
        return False

    # Compare non-NaN values
    valid_mask = ~ref_nan_mask
    ref_valid = ref_data[valid_mask]
    test_valid = test_data[valid_mask]

    if not np.allclose(ref_valid, test_valid, rtol=tolerance, atol=tolerance):
        diff = np.abs(ref_valid - test_valid)
        max_diff = diff.max()
        mean_diff = diff.mean()

        logger.error(f"Data mismatch for '{name}':")
        logger.error(f"  Max difference: {max_diff}")
        logger.error(f"  Mean difference: {mean_diff}")
        logger.error(f"  Tolerance: {tolerance}")
        return False

    logger.info(f"  ✓ {name}")
    return True
