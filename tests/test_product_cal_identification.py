"""DISP-CAL product: identification metadata taken from the DISP input."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pytest
import xarray as xr

from cal_disp.product import CalProduct, DispProduct
from cal_disp.product._disp import read_disp_metadata
from cal_disp.product.output._identification import iso_utc
from cal_disp.product.output._utils import (
    bounding_polygon_wkt,
    geotransform_from_coords,
    grid_bounds,
)
from cal_disp.workflow import run_calibration

# /identification values as in the real DISP-S1 F08882 product
DISP_IDENTIFICATION = {
    "absolute_orbit_number": 41406,
    "track_number": 34,
    "orbit_pass_direction": "Ascending",
    "look_direction": "Right",
    "instrument_name": "C-SAR",
    "radar_band": "C",
    "acquisition_mode": "IW",
    "source_data_satellite_names": "S1A",
    "source_data_dem_name": "Copernicus GLO-30",
    "source_data_imaging_geometry": "Geocoded",
    "product_data_access": (
        "https://search.asf.alaska.edu/#/?dataset=OPERA-S1&productTypes=DISP-S1"
    ),
    "static_layers_data_access": (
        "https://search.asf.alaska.edu/#/?dataset=OPERA-S1"
        "&productTypes=DISP-S1-STATIC&frame=8882"
    ),
    "ceos_analysis_ready_data_document_identifier": (
        "https://ceos.org/ard/files/PFS/SAR/v1.2/"
        "CEOS-ARD_PFS_Synthetic_Aperture_Radar_v1.2.pdf"
    ),
    "ceos_analysis_ready_data_product_type": "InSAR",
}
# Deliberately not the old hard-coded "S1A"
DISP_METADATA = {"platform_id": "S1B"}


@pytest.fixture
def disp_with_identification(sample_disp_product: Path) -> Path:
    """The synthetic DISP product with a DISP-S1-like identification group."""
    with h5py.File(sample_disp_product, "a") as f:
        ident = f["identification"]
        for name, value in DISP_IDENTIFICATION.items():
            ident.create_dataset(name, data=value)
        meta = f.require_group("metadata")
        for name, value in DISP_METADATA.items():
            meta.create_dataset(name, data=value)
    return sample_disp_product


def _scalars(ds: xr.Dataset) -> dict:
    return {str(k): v.item() for k, v in ds.items()}


def test_read_disp_metadata(disp_with_identification: Path, caplog):
    meta = read_disp_metadata(disp_with_identification)

    assert meta["track_number"] == 34
    assert meta["absolute_orbit_number"] == 41406
    assert meta["orbit_pass_direction"] == "Ascending"
    assert meta["platform_id"] == "S1B"
    assert meta["radar_wavelength"] == pytest.approx(0.05546)
    # Missing in the synthetic product: None, with a warning naming the path
    assert meta["bounding_polygon"] is None
    assert "/identification/bounding_polygon" in caplog.text
    assert DispProduct.from_path(disp_with_identification).read_metadata() == meta


def test_identification_from_disp(
    tmp_path: Path,
    disp_with_identification: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    lookup_file, tenv8_dir = sample_unr_data
    out = run_calibration(
        disp_file=disp_with_identification,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
    )
    product = CalProduct.from_path(out)
    ident = _scalars(product.open_identification())
    meta = _scalars(product.open_metadata())

    assert ident["track_number"] == 34
    assert ident["absolute_orbit_number"] == 41406
    assert ident["orbit_pass_direction"] == "ascending"
    assert ident["look_direction"] == "right"
    assert ident["instrument_name"] == "C-SAR"
    assert ident["radar_band"] == "C"
    assert ident["source_data_satellite_names"] == "S1A"
    assert ident["source_data_dem_name"] == "Copernicus GLO-30"
    assert ident["source_data_imaging_geometry"] == "Geocoded"
    assert ident["static_layers_data_access"] == (
        DISP_IDENTIFICATION["static_layers_data_access"]
    )
    # The source of a DISP-CAL product is the DISP product
    assert ident["source_data_access"] == DISP_IDENTIFICATION["product_data_access"]
    assert ident["ceos_analysis_ready_data_product_type"] == "InSAR"
    assert ident["ceos_analysis_ready_data_document_identifier"].startswith(
        "https://ceos.org"
    )
    assert "asf.alaska.edu" in ident["product_data_access"]
    assert ident["processing_facility"] == "NASA Jet Propulsion Laboratory on AWS"
    assert meta["platform_id"] == "S1B"
    for value in ident.values():
        assert "example.com" not in str(value)

    # Datetimes: ISO 8601 UTC with Z
    assert ident["reference_datetime"] == "2022-01-11T00:26:51Z"
    assert ident["secondary_datetime"] == "2022-07-22T00:26:57Z"
    datetime.strptime(ident["processing_start_datetime"], "%Y-%m-%dT%H:%M:%SZ")

    # Geometry: outer pixel edges (UTM) and a lon/lat polygon
    assert ident["product_bounding_box"] == "(404950.0, 3758050.0, 424950.0, 3778050.0)"
    poly = ident["bounding_polygon"]
    assert poly.startswith("POLYGON((") and poly.endswith("))")
    lons, lats = zip(
        *(map(float, p.split()) for p in poly[len("POLYGON((") : -2].split(", "))
    )
    assert -118.1 < min(lons) < max(lons) < -117.7
    assert 33.9 < min(lats) < max(lats) < 34.2
    assert ident["product_sample_spacing"] == "100m"

    # nodata: the product's own NaNs
    with xr.open_dataset(out) as ds:
        assert ident["nodata_pixel_count"] == int(np.isnan(ds["calibration"]).sum())
    # every staged station file of the grid type is listed (no truncation)
    listed = ident["source_calibration_file_list"].split(", ")
    staged = {p.name for p in tenv8_dir.glob("*_IGS20_constant.tenv8")}
    assert listed[0] == lookup_file.name
    assert len(staged) > 1
    assert set(listed[1:]) == staged


def test_access_urls_from_arguments_override_disp(
    tmp_path: Path,
    disp_with_identification: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    lookup_file, tenv8_dir = sample_unr_data
    out = run_calibration(
        disp_file=disp_with_identification,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        processing_facility="Test facility",
        product_data_access="https://doi.org/10.0000/cal",
        static_layers_data_access="https://doi.org/10.0000/static",
        source_data_access="https://doi.org/10.0000/disp",
    )
    ident = _scalars(CalProduct.from_path(out).open_identification())

    assert ident["processing_facility"] == "Test facility"
    assert ident["product_data_access"] == "https://doi.org/10.0000/cal"
    assert ident["static_layers_data_access"] == "https://doi.org/10.0000/static"
    assert ident["source_data_access"] == "https://doi.org/10.0000/disp"


def test_identification_fallbacks_when_disp_lacks_fields(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    caplog,
):
    """No stubbed S1A/ascending/orbit 0: marked fallbacks and a WARNING."""
    lookup_file, tenv8_dir = sample_unr_data
    out = run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
    )
    product = CalProduct.from_path(out)
    ident = _scalars(product.open_identification())

    assert ident["track_number"] == -1
    assert ident["absolute_orbit_number"] == -1
    assert ident["orbit_pass_direction"] == "unknown"
    assert ident["source_data_satellite_names"] == "unknown"
    assert ident["static_layers_data_access"] == "unknown"
    assert ident["source_data_access"] == "unknown"
    assert "ceos_analysis_ready_data_product_type" not in ident
    assert _scalars(product.open_metadata())["platform_id"] == "unknown"
    assert "/identification/track_number" in caplog.text
    assert any(
        r.levelname == "WARNING" and r.name == "cal_disp.product._disp"
        for r in caplog.records
    )


def test_pge_runconfig_embedded_verbatim(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    sample_algorithm_params: Path,
):
    """/metadata/pge_runconfig is the RunConfig YAML that was run."""
    from ruamel.yaml import YAML

    from cal_disp.cli.config import create_config
    from cal_disp.cli.run import run_main

    lookup_file, tenv8_dir = sample_unr_data
    config_file = create_config(
        disp_file=sample_disp_product,
        frame_id=8882,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        unr_grid_version="0.2",
        unr_grid_type="constant",
        algorithm_params_file=sample_algorithm_params,
        los_file=sample_static_los,
        dem_file=sample_static_dem,
        output_dir=tmp_path / "out",
        work_dir=tmp_path / "work",
    )
    out = run_main(config_file)

    embedded = _scalars(CalProduct.from_path(out).open_metadata())["pge_runconfig"]
    assert embedded == config_file.read_text()
    assert not embedded.startswith("{")  # not a Python dict repr
    parsed = YAML(typ="safe").load(embedded)["cal_disp_workflow"]
    assert parsed["input_file_group"]["frame_id"] == 8882
    assert parsed["output_options"]["compression"] is True
    assert "product_data_access" in parsed["output_options"]


def test_iso_utc():
    assert iso_utc(datetime(2022, 1, 11, 0, 26, 51)) == "2022-01-11T00:26:51Z"
    aware = datetime(2022, 1, 11, 1, 26, 51, 123, tzinfo=timezone.utc)
    assert iso_utc(aware) == "2022-01-11T01:26:51Z"


def test_grid_geometry_helpers():
    # The real F08882 grid: centres start at 71985, DISP GeoTransform at 71970
    x = 71985.0 + 30.0 * np.arange(4)
    y = 3385905.0 - 30.0 * np.arange(3)
    assert geotransform_from_coords(x, y) == "71970.0 30.0 0.0 3385920.0 0.0 -30.0"
    assert grid_bounds(x, y) == (71970.0, 3385830.0, 72090.0, 3385920.0)
    wkt = bounding_polygon_wkt(x, y, "EPSG:32615", points_per_edge=2)
    assert wkt.startswith("POLYGON((") and wkt.count(",") == 8
    lon, lat = map(float, wkt[len("POLYGON((") :].split(",")[0].split())
    assert -98 < lon < -96 and 30 < lat < 31
