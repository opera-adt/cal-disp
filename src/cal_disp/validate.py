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
        Which group to validate: "main", "model_3d", or "all". Default is "all".

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

    # Compare main group
    if group in ("main", "all"):
        if not _compare_group(ref_cal, test_cal, "main", tolerance):
            return False

    # Compare model_3d group if it exists
    if group in ("model_3d", "all"):
        ref_has_model = ref_cal.has_model_3d()
        test_has_model = test_cal.has_model_3d()

        if ref_has_model != test_has_model:
            logger.error(
                f"model_3d group mismatch: reference={ref_has_model},"
                f" test={test_has_model}"
            )
            return False

        if ref_has_model and not _compare_group(
            ref_cal, test_cal, "model_3d", tolerance
        ):
            return False

    logger.info("✓ Validation passed")
    return True


def _compare_metadata(ref: CalProduct, test: CalProduct) -> bool:
    """Compare product metadata."""
    checks = {
        "frame_id": (ref.frame_id, test.frame_id),
        "sensor": (ref.sensor, test.sensor),
        "mode": (ref.mode, test.mode),
        "polarization": (ref.polarization, test.polarization),
        "primary_date": (ref.primary_date, test.primary_date),
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
    group_name = "main" if group == "main" else "model_3d"
    logger.info(f"Validating {group_name} group...")

    # Open datasets
    if group == "main":
        ref_ds = ref.open_dataset()
        test_ds = test.open_dataset()
    else:
        ref_ds = ref.open_model_3d()
        test_ds = test.open_model_3d()

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
