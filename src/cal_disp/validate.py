"""Compare a DISP-CAL product against a reference (golden) product.

Everything that makes a product is compared, group by group:

* presence of every group (``main``, ``identification``, ``metadata``,
  ``auxiliary``); a group missing on one side is a failure,
* the set of variables and coordinates, and their on-disk dtype and shape,
* variable and group attributes (``units``, ``long_name``, ``grid_mapping``,
  ``_FillValue``, ...), except the volatile fields in :data:`IGNORED_ATTRS`,
* the CRS (``spatial_ref``, compared as a :class:`rasterio.crs.CRS`) and the
  grid transform (``GeoTransform`` attribute and the ``x``/``y`` coordinates),
* values: strings and integers exactly, floating point with ``rtol``/``atol``
  (reference zeros included: the NaN pattern must match and every finite
  value is compared),
* the browse image next to the test product (must be at most
  :data:`BROWSE_MAX_DIM` pixels on each side; only a warning when absent).

Differences are collected and logged grouped by section before the result is
returned, so one run shows everything that is wrong.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    import xarray as xr
    from rasterio.crs import CRS

    from cal_disp.product import CalProduct

logger = logging.getLogger(__name__)

#: Default relative and absolute tolerance for floating point comparison.
#: The CLI, the scripts and :func:`compare_cal_products` all use this value.
DEFAULT_TOLERANCE: float = 1e-6

#: Largest allowed browse image side (pixels).
BROWSE_MAX_DIM: int = 2048

#: Groups of a DISP-CAL product; ``main`` is the root group.
PRODUCT_GROUPS: tuple[str, ...] = ("main", "identification", "metadata", "auxiliary")

#: Attributes and scalar variables that legitimately differ between two runs
#: of the same software on the same inputs. Keys are ``(group, variable)``;
#: ``variable`` is ``""`` for group attributes and ``"*"`` for every variable
#: of the group. Keep this list short: everything not listed is compared.
IGNORED_ATTRS: dict[tuple[str, str], frozenset[str]] = {
    # Package version of the build that wrote the file, the HDF5/h5netcdf
    # library versions recorded by the writer, and the creation-time
    # ``history`` line (timestamp + version).
    ("main", ""): frozenset({"software_version", "_NCProperties", "history"}),
}

#: Scalar variables whose *value* is volatile: production time, build
#: versions, and the embedded runconfig (it carries the absolute input and
#: scratch paths of the run). Their attributes are still compared.
IGNORED_VALUES: dict[str, frozenset[str]] = {
    "identification": frozenset({"processing_start_datetime"}),
    "metadata": frozenset(
        {
            "cal_disp_software_version",
            "venti_software_version",
            "pge_runconfig",
        }
    ),
}

# Attributes of the CRS variable that are compared as a CRS, not as text
_CRS_TEXT_ATTRS = frozenset({"crs_wkt", "spatial_ref"})
# netCDF-4 bookkeeping written by the HDF5 layer, not part of the product
_HDF5_INTERNAL_ATTRS = frozenset(
    {"DIMENSION_LIST", "REFERENCE_LIST", "CLASS", "NAME", "_Netcdf4Dimid"}
)


@dataclass
class ValidationReport:
    """Differences found while comparing two products, grouped by section."""

    failures: dict[str, list[str]] = field(default_factory=dict)
    warnings: dict[str, list[str]] = field(default_factory=dict)

    def fail(self, section: str, message: str) -> None:
        """Record a difference (the validation fails)."""
        self.failures.setdefault(section, []).append(message)

    def warn(self, section: str, message: str) -> None:
        """Record a warning (the validation still passes)."""
        self.warnings.setdefault(section, []).append(message)

    @property
    def passed(self) -> bool:
        """True when no difference was recorded."""
        return not self.failures

    @property
    def n_failures(self) -> int:
        """Number of recorded differences."""
        return sum(len(v) for v in self.failures.values())

    def log(self) -> None:
        """Log every warning and failure, grouped by section."""
        for section, messages in self.warnings.items():
            logger.warning("[%s]", section)
            for msg in messages:
                logger.warning("  %s", msg)
        for section, messages in self.failures.items():
            logger.error("[%s] %d difference(s)", section, len(messages))
            for msg in messages:
                logger.error("  %s", msg)


def compare_cal_products(
    reference_file: Path,
    test_file: Path,
    rtol: float = DEFAULT_TOLERANCE,
    atol: float = DEFAULT_TOLERANCE,
    group: str = "all",
    tolerance: float | None = None,
    browse_max_dim: int = BROWSE_MAX_DIM,
) -> bool:
    """Compare two DISP-CAL products.

    Parameters
    ----------
    reference_file : Path
        Reference (golden) product file.
    test_file : Path
        Test product file to validate.
    rtol, atol : float, optional
        Relative and absolute tolerance for floating point values
        (``numpy.allclose`` semantics). Default :data:`DEFAULT_TOLERANCE`.
    group : str, optional
        ``"main"`` (root, identification and metadata groups and the browse
        image), ``"auxiliary"`` (the 3-D model group) or ``"all"``.
    tolerance : float, optional
        Shortcut setting both ``rtol`` and ``atol``.
    browse_max_dim : int, optional
        Largest allowed browse image side in pixels.

    Returns
    -------
    bool
        True if the products match, False otherwise. All differences are
        logged before returning.

    """
    from cal_disp.product import CalProduct

    if tolerance is not None:
        rtol = atol = tolerance
    if group not in ("main", "auxiliary", "all"):
        raise ValueError(f"group must be 'main', 'auxiliary' or 'all', got {group!r}")

    logger.info(
        "Comparing %s against %s (rtol=%g, atol=%g, group=%s)",
        test_file.name,
        reference_file.name,
        rtol,
        atol,
        group,
    )
    report = ValidationReport()

    ref_cal = CalProduct.from_path(reference_file)
    test_cal = CalProduct.from_path(test_file)
    _compare_filename_metadata(ref_cal, test_cal, report)

    groups = list(PRODUCT_GROUPS)
    if group == "main":
        groups.remove("auxiliary")
    elif group == "auxiliary":
        groups = ["auxiliary"]

    for name in groups:
        ref_ds = _open_group(reference_file, name)
        test_ds = _open_group(test_file, name)
        if ref_ds is None and test_ds is None:
            if name == "auxiliary":
                report.warn(name, "group absent from both products")
            else:
                report.fail(name, "group absent from both products")
            continue
        if ref_ds is None or test_ds is None:
            missing = "test" if test_ds is None else "reference"
            report.fail(name, f"group missing from the {missing} product")
            continue
        with ref_ds, test_ds:
            if name == "main":
                _validate_product_structure(test_ds, report)
            if name == "identification":
                _validate_identification_structure(test_ds, report)
            _compare_group(ref_ds, test_ds, name, rtol, atol, report)

    if group != "auxiliary":
        _check_browse_image(test_file, browse_max_dim, report)

    report.log()
    if report.passed:
        logger.info("Validation passed")
    else:
        logger.error("Validation FAILED: %d difference(s)", report.n_failures)
    return report.passed


# Group access


def _open_group(path: Path, name: str) -> xr.Dataset | None:
    """Open one product group raw (no CF decoding), or None if it is absent."""
    import xarray as xr

    kwargs: dict[str, Any] = {"engine": "h5netcdf", "decode_cf": False}
    if name != "main":
        kwargs["group"] = name
    try:
        return xr.open_dataset(path, **kwargs)
    except (OSError, ValueError, KeyError):
        return None


# Filename and structure checks


def _compare_filename_metadata(
    ref: CalProduct, test: CalProduct, report: ValidationReport
) -> None:
    """Compare the fields encoded in the product filename (not production time)."""
    checks = {
        "frame_id": (ref.frame_id, test.frame_id),
        "sensor": (ref.sensor, test.sensor),
        "mode": (ref.mode, test.mode),
        "polarization": (ref.polarization, test.polarization),
        "reference_date": (ref.reference_date, test.reference_date),
        "secondary_date": (ref.secondary_date, test.secondary_date),
        "version": (ref.version, test.version),
    }
    for name, (ref_val, test_val) in checks.items():
        if ref_val != test_val:
            report.fail("filename", f"{name}: reference={ref_val}, test={test_val}")


def _validate_product_structure(ds: xr.Dataset, report: ValidationReport) -> None:
    """Check dims (calibration is (time, y, x), none unnamed)."""
    if "calibration" in ds.data_vars:
        cal_dims = ds["calibration"].dims
        if cal_dims != ("time", "y", "x"):
            report.fail(
                "main", f"calibration has dimensions {cal_dims}, expected (time, y, x)"
            )
    for var in ds.data_vars:
        if any("dim_" in str(dim) for dim in ds[var].dims):
            report.fail("main", f"'{var}' has unnamed dimension: {ds[var].dims}")


def _validate_identification_structure(
    ds: xr.Dataset, report: ValidationReport
) -> None:
    """Check that the identification list variables are scalar strings."""
    for var in (
        "source_calibration_file_list",
        "source_data_file_list",
        "source_data_satellite_names",
    ):
        if var not in ds.data_vars:
            continue
        data_arr = ds[var]
        if data_arr.dims != ():
            report.fail(
                "identification",
                f"'{var}' should be scalar but has dimensions {data_arr.dims}",
            )
            continue
        value = data_arr.item()
        if isinstance(value, bytes):
            value = value.decode()
        if not isinstance(value, str):
            report.fail(
                "identification", f"'{var}' should be a string, got {type(value)}"
            )


# Group comparison


def _compare_group(
    ref_ds: xr.Dataset,
    test_ds: xr.Dataset,
    group: str,
    rtol: float,
    atol: float,
    report: ValidationReport,
) -> None:
    """Compare attributes, variables, CRS and transform of one group."""
    _compare_attrs(
        ref_ds.attrs,
        test_ds.attrs,
        group,
        "group attribute",
        IGNORED_ATTRS.get((group, ""), frozenset()),
        rtol,
        atol,
        report,
    )

    ref_vars = set(ref_ds.variables)
    test_vars = set(test_ds.variables)
    if ref_vars != test_vars:
        if ref_vars - test_vars:
            report.fail(group, f"missing from test: {sorted(ref_vars - test_vars)}")
        if test_vars - ref_vars:
            report.fail(group, f"only in test: {sorted(test_vars - ref_vars)}")

    for var in sorted(ref_vars & test_vars):
        if var == "spatial_ref":
            _compare_spatial_ref(ref_ds, test_ds, group, rtol, atol, report)
            continue
        ignored = IGNORED_ATTRS.get((group, "*"), frozenset()) | IGNORED_ATTRS.get(
            (group, var), frozenset()
        )
        _compare_variable(
            ref_ds[var],
            test_ds[var],
            group,
            var,
            rtol,
            atol,
            report,
            ignored_attrs=ignored,
            compare_values=var not in IGNORED_VALUES.get(group, frozenset()),
        )

    if "x" in ref_ds.coords and "y" in ref_ds.coords:
        _compare_grid_from_coords(ref_ds, test_ds, group, report)


def _compare_variable(
    ref: xr.DataArray,
    test: xr.DataArray,
    group: str,
    name: str,
    rtol: float,
    atol: float,
    report: ValidationReport,
    ignored_attrs: frozenset[str] = frozenset(),
    compare_values: bool = True,
) -> None:
    """Compare dims, dtype, shape, attributes and values of one variable."""
    label = f"{group}/{name}"
    if ref.dims != test.dims:
        report.fail(group, f"{name}: dims {test.dims}, expected {ref.dims}")
        return
    if ref.shape != test.shape:
        report.fail(group, f"{name}: shape {test.shape}, expected {ref.shape}")
        return

    ref_dtype, test_dtype = _disk_dtype(ref), _disk_dtype(test)
    if ref_dtype != test_dtype:
        report.fail(group, f"{name}: dtype {test_dtype}, expected {ref_dtype}")

    _compare_attrs(
        ref.attrs,
        test.attrs,
        group,
        f"{name} attribute",
        ignored_attrs,
        rtol,
        atol,
        report,
    )

    if not compare_values:
        logger.debug("Skipping volatile value of %s", label)
        return
    _compare_values(ref.values, test.values, group, name, rtol, atol, report)


def _disk_dtype(da: xr.DataArray) -> str:
    """On-disk dtype of a variable ('str' for variable-length strings)."""
    dtype = da.encoding.get("dtype", da.dtype)
    if np.dtype(dtype).kind in ("O", "U", "S"):
        return "str"
    return str(np.dtype(dtype))


def _compare_values(
    ref_data: np.ndarray,
    test_data: np.ndarray,
    group: str,
    name: str,
    rtol: float,
    atol: float,
    report: ValidationReport,
) -> None:
    """Compare array values: exact for strings/integers, tolerance for floats."""
    ref_data = np.asarray(ref_data)
    test_data = np.asarray(test_data)

    if ref_data.dtype.kind in ("O", "U", "S") or test_data.dtype.kind in (
        "O",
        "U",
        "S",
    ):
        ref_s, test_s = _as_str_array(ref_data), _as_str_array(test_data)
        if not np.array_equal(ref_s, test_s):
            report.fail(
                group,
                f"{name}: {_short(test_s)!r}, expected {_short(ref_s)!r}",
            )
        return

    if ref_data.dtype.kind in ("i", "u", "b"):
        if not np.array_equal(ref_data, test_data):
            report.fail(
                group,
                (
                    f"{name}: {_short(test_data)}, expected {_short(ref_data)}"
                    if ref_data.ndim == 0
                    else f"{name}: {int((ref_data != test_data).sum())} value(s) differ"
                ),
            )
        return

    ref_nan = np.isnan(ref_data)
    test_nan = np.isnan(test_data)
    if not np.array_equal(ref_nan, test_nan):
        report.fail(
            group,
            f"{name}: NaN pattern differs (reference {int(ref_nan.sum())} NaN,"
            f" test {int(test_nan.sum())} NaN)",
        )
        return

    # Every finite value is compared, reference zeros included
    valid = ~ref_nan
    ref_valid = ref_data[valid]
    test_valid = test_data[valid]
    close = np.isclose(ref_valid, test_valid, rtol=rtol, atol=atol)
    if close.all():
        return
    diff = np.abs(ref_valid.astype(np.float64) - test_valid.astype(np.float64))
    n_bad = int((~close).sum())
    if ref_data.ndim == 0:
        report.fail(
            group,
            f"{name}: {test_valid[0]!r}, expected {ref_valid[0]!r}"
            f" (|diff|={diff[0]:.3g}, rtol={rtol:g}, atol={atol:g})",
        )
        return
    first = np.argwhere(~close)[0][0]
    first_idx = tuple(int(i) for i in np.argwhere(valid)[first])
    report.fail(
        group,
        f"{name}: {n_bad} of {valid.sum()} values beyond tolerance"
        f" (max |diff|={diff.max():.3g}, mean |diff|={diff.mean():.3g},"
        f" rtol={rtol:g}, atol={atol:g}); first at index {first_idx}:"
        f" test={test_valid[first]!r}, reference={ref_valid[first]!r}",
    )


def _as_str_array(data: np.ndarray) -> np.ndarray:
    """Decode bytes and return a unicode string array."""
    flat = [v.decode() if isinstance(v, bytes) else str(v) for v in data.ravel()]
    return np.array(flat, dtype=str).reshape(data.shape)


def _short(data: np.ndarray, limit: int = 80) -> Any:
    """Scalar value for messages, shortened if it is a long string."""
    value = data.item() if data.ndim == 0 else data
    if isinstance(value, str) and len(value) > limit:
        return value[: limit - 3] + "..."
    return value


# Attributes


def _compare_attrs(
    ref_attrs: dict[str, Any],
    test_attrs: dict[str, Any],
    group: str,
    what: str,
    ignored: frozenset[str],
    rtol: float,
    atol: float,
    report: ValidationReport,
) -> None:
    """Compare two attribute dicts (numbers with tolerance, everything else exactly)."""
    ref_keys = set(ref_attrs) - ignored - _HDF5_INTERNAL_ATTRS
    test_keys = set(test_attrs) - ignored - _HDF5_INTERNAL_ATTRS
    for key in sorted(ref_keys - test_keys):
        report.fail(group, f"{what} '{key}' missing from test")
    for key in sorted(test_keys - ref_keys):
        report.fail(group, f"{what} '{key}' only in test: {test_attrs[key]!r}")
    for key in sorted(ref_keys & test_keys):
        if not _attr_equal(ref_attrs[key], test_attrs[key], rtol, atol):
            report.fail(
                group,
                f"{what} '{key}': {_short_attr(test_attrs[key])!r},"
                f" expected {_short_attr(ref_attrs[key])!r}",
            )


def _attr_equal(a: Any, b: Any, rtol: float, atol: float) -> bool:
    a_arr, b_arr = np.asarray(a), np.asarray(b)
    if a_arr.dtype.kind in "fiub" and b_arr.dtype.kind in "fiub":
        if a_arr.shape != b_arr.shape:
            return False
        return bool(np.allclose(a_arr, b_arr, rtol=rtol, atol=atol, equal_nan=True))
    if a_arr.dtype.kind in "fiub" or b_arr.dtype.kind in "fiub":
        return False
    return bool(np.array_equal(_as_str_array(a_arr), _as_str_array(b_arr)))


def _short_attr(value: Any, limit: int = 80) -> Any:
    if isinstance(value, bytes):
        value = value.decode(errors="replace")
    if isinstance(value, str) and len(value) > limit:
        return value[: limit - 3] + "..."
    return value


# CRS and transform


def _compare_spatial_ref(
    ref_ds: xr.Dataset,
    test_ds: xr.Dataset,
    group: str,
    rtol: float,
    atol: float,
    report: ValidationReport,
) -> None:
    """Compare the CRS (as rasterio CRS) and the GeoTransform of a group."""
    from rasterio.crs import CRS

    ref_attrs, test_attrs = ref_ds["spatial_ref"].attrs, test_ds["spatial_ref"].attrs

    ref_wkt = ref_attrs.get("crs_wkt", ref_attrs.get("spatial_ref"))
    test_wkt = test_attrs.get("crs_wkt", test_attrs.get("spatial_ref"))
    if ref_wkt is None or test_wkt is None:
        missing = "test" if test_wkt is None else "reference"
        report.fail(group, f"spatial_ref: no crs_wkt in the {missing} product")
    else:
        try:
            ref_crs, test_crs = CRS.from_wkt(ref_wkt), CRS.from_wkt(test_wkt)
        except Exception as exc:  # rasterio raises CRSError
            report.fail(group, f"spatial_ref: unreadable crs_wkt ({exc})")
        else:
            if ref_crs != test_crs:
                report.fail(
                    group,
                    f"spatial_ref: CRS {_crs_name(test_crs)},"
                    f" expected {_crs_name(ref_crs)}",
                )

    ref_gt, test_gt = ref_attrs.get("GeoTransform"), test_attrs.get("GeoTransform")
    if (ref_gt is None) != (test_gt is None):
        missing = "test" if test_gt is None else "reference"
        report.fail(group, f"spatial_ref: GeoTransform missing from the {missing}")
    elif ref_gt is not None:
        ref_vals = np.array([float(v) for v in str(ref_gt).split()])
        test_vals = np.array([float(v) for v in str(test_gt).split()])
        if ref_vals.shape != test_vals.shape or not np.allclose(
            ref_vals, test_vals, rtol=0, atol=_transform_atol(ref_vals[1])
        ):
            report.fail(
                group, f"spatial_ref: GeoTransform '{test_gt}', expected '{ref_gt}'"
            )

    # Remaining attributes (ellipsoid, projection parameters, ...) as usual
    _compare_attrs(
        {k: v for k, v in ref_attrs.items() if k not in _CRS_TEXT_ATTRS},
        {k: v for k, v in test_attrs.items() if k not in _CRS_TEXT_ATTRS},
        group,
        "spatial_ref attribute",
        frozenset({"GeoTransform"}),
        rtol,
        atol,
        report,
    )


def _crs_name(crs: CRS) -> str:
    epsg = crs.to_epsg()
    return f"EPSG:{epsg}" if epsg else crs.to_string()[:60]


def _transform_atol(pixel_size: float) -> float:
    """Tolerance for transform terms: 1e-6 of a pixel (at least 1e-9)."""
    return max(abs(float(pixel_size)) * 1e-6, 1e-9)


def _compare_grid_from_coords(
    ref_ds: xr.Dataset, test_ds: xr.Dataset, group: str, report: ValidationReport
) -> None:
    """Compare origin and spacing derived from the x/y coordinates."""
    for coord in ("x", "y"):
        ref_c = np.asarray(ref_ds[coord].values, dtype=np.float64)
        test_c = np.asarray(test_ds[coord].values, dtype=np.float64)
        if ref_c.shape != test_c.shape or ref_c.size < 2:
            continue  # shape difference already reported
        ref_step, test_step = ref_c[1] - ref_c[0], test_c[1] - test_c[0]
        tol = _transform_atol(ref_step)
        if not (
            math.isclose(ref_c[0], test_c[0], rel_tol=0, abs_tol=tol)
            and math.isclose(ref_step, test_step, rel_tol=0, abs_tol=tol)
        ):
            report.fail(
                group,
                f"{coord} grid: origin {test_c[0]}, step {test_step}; expected"
                f" origin {ref_c[0]}, step {ref_step}",
            )


# Browse image


def _check_browse_image(
    test_file: Path, max_dim: int, report: ValidationReport
) -> None:
    """Check the browse PNG next to the product: present and at most max_dim px."""
    png = test_file.with_suffix(".png")
    if not png.exists():
        report.warn("browse", f"no browse image next to the test product: {png.name}")
        return
    try:
        from PIL import Image

        with Image.open(png) as img:
            width, height = img.size
    except Exception as exc:
        report.fail("browse", f"cannot read browse image {png.name}: {exc}")
        return
    if max(width, height) > max_dim:
        report.fail(
            "browse",
            f"{png.name} is {width}x{height} px; each side must be <= {max_dim} px",
        )
    else:
        logger.info("Browse image %s: %dx%d px", png.name, width, height)
