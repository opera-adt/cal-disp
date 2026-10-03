"""DISP-CAL product: ``spatial_ref`` keeps its GeoTransform; GeoTIFF export."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
import xarray as xr
from cal_product_helpers import write_product

from cal_disp.product import CalProduct
from cal_disp.product.output._utils import get_transform, make_spatial_ref


def test_make_spatial_ref_keeps_disp_attrs():
    disp_attrs = {
        "crs_wkt": rasterio.crs.CRS.from_epsg(32615).to_wkt(),
        "grid_mapping_name": "transverse_mercator",
        "semi_major_axis": np.array([6378137.0]),
        "GeoTransform": "71970.0 30.0 0.0 3385920.0 0.0 -30.0",
        "units": "unitless",
        "long_name": "Dummy variable with geo-referencing metadata in attributes",
    }
    x = 71985.0 + 30.0 * np.arange(10)
    y = 3385905.0 - 30.0 * np.arange(8)

    out = make_spatial_ref(xr.DataArray(0, attrs=disp_attrs), x, y)

    # Same grid as DISP: every attribute is preserved
    assert out.attrs["GeoTransform"] == disp_attrs["GeoTransform"]
    assert out.attrs["crs_wkt"] == disp_attrs["crs_wkt"]
    assert out.attrs["units"] == "unitless"
    assert out.attrs["long_name"] == disp_attrs["long_name"]
    assert out.attrs["semi_major_axis"] == 6378137.0  # scalar, not a 1-array
    # A coarser grid gets its own GeoTransform
    coarse = make_spatial_ref(xr.DataArray(0, attrs=disp_attrs), x[::5], y[::4])
    assert coarse.attrs["GeoTransform"] == "71910.0 150.0 0.0 3385965.0 0.0 -120.0"


def test_spatial_ref_and_geotiff(tmp_path: Path, sample_disp_product: Path):
    product, cal = write_product(sample_disp_product, tmp_path / "out")
    with xr.open_dataset(sample_disp_product, engine="h5netcdf") as disp:
        disp_attrs = dict(disp["spatial_ref"].attrs)

    # GeoTransform follows the DISP-S1 convention: x/y are pixel centres, the
    # transform origin is the outer edge (half a pixel out). The synthetic
    # fixture's own GeoTransform puts the origin on the first centre instead;
    # the real product (F08882: x[0]=71985, origin 71970) is consistent.
    for group, expected in (
        (None, "404950.0 100.0 0.0 3778050.0 0.0 -100.0"),
        # coarse grid (every 50th pixel): 5 km pixels
        ("auxiliary", "402500.0 5000.0 0.0 3780500.0 0.0 -5000.0"),
    ):
        with xr.open_dataset(product.path, group=group, engine="h5netcdf") as ds:
            attrs = ds["spatial_ref"].attrs
            assert attrs["crs_wkt"] == disp_attrs["crs_wkt"]
            assert attrs["GeoTransform"] == expected
            assert attrs["units"] == "unitless"
            assert attrs["grid_mapping_name"] == "transverse_mercator"
            assert get_transform(ds) == rasterio.transform.from_origin(
                *(float(expected.split()[i]) for i in (0, 3, 1)),
                -float(expected.split()[5]),
            )

    tif = CalProduct.from_path(product.path).to_geotiff(
        "calibration", tmp_path / "cal.tif"
    )
    with rasterio.open(tif) as src:
        assert src.crs.to_epsg() == 32611
        assert src.transform == rasterio.transform.from_origin(
            404950.0, 3778050.0, 100.0, 100.0
        )
        assert src.read(1).tobytes() == cal[0].tobytes()
        assert src.tags()["frame_id"] == "8882"

    aux_tif = CalProduct.from_path(product.path).to_geotiff(
        "up_down", tmp_path / "up.tif", group="auxiliary", tiled=False
    )
    with rasterio.open(aux_tif) as src:
        assert src.transform.a == 5000.0 and src.shape == (4, 4)
