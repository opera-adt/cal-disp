"""End-to-end test of ``run_calibration`` against the installed Venti."""

from __future__ import annotations

import socket
from pathlib import Path

import numpy as np
import pytest
import rasterio
import xarray as xr
from _gamma import gamma_parameters

from cal_disp.config._algorithm import AlgorithmParameters, CalibrationOptions
from cal_disp.workflow import run_calibration


@pytest.fixture
def no_network(monkeypatch):
    """Fail on any network connection (GNSS must come from pre-staged files)."""

    def _blocked(*_args, **_kwargs):
        raise RuntimeError("network access attempted during an offline run")

    monkeypatch.setattr(socket.socket, "connect", _blocked)
    monkeypatch.setattr(socket, "create_connection", _blocked)


@pytest.mark.usefixtures("no_network")
@pytest.mark.parametrize("grid_type", ["constant", "variable"])
def test_run_calibration_offline(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    grid_type: str,
):
    lookup_file, tenv8_dir = sample_unr_data
    output_dir = tmp_path / "out"

    out_path = run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=output_dir,
        los_file=sample_static_los,
        algorithm_parameters=gamma_parameters(
            {"calibration_options": {"grid_type": grid_type}}
        ),
    )

    assert out_path.exists()
    with xr.open_dataset(out_path) as ds:
        calibration = ds["calibration"].values
        calibration_std = ds["calibration_std"].values
        uncertainty_source = ds["calibration_std"].attrs["uncertainty_source"]
    with xr.open_dataset(sample_disp_product) as disp:
        valid = np.isfinite(disp["displacement"].values)

    assert calibration.shape == (1, *valid.shape)
    assert np.isfinite(calibration[0][valid]).all()
    # calibration_std holds the GNSS reference uncertainty (m), documented as such
    assert calibration_std.shape == calibration.shape
    assert (calibration_std[0][valid] > 0).all()
    assert (calibration_std[0][valid] < 0.01).all()
    assert uncertainty_source.startswith("GNSS reference uncertainty")
    # GNSS came only from the pre-staged files of this grid type
    staged = {p.name for p in tenv8_dir.glob(f"*_IGS20_{grid_type}.tenv8")}
    used = {p.name for p in (output_dir / "scratch" / "gnss").glob("*.tenv8")}
    assert used == staged
    assert len(staged) == 4


def test_run_calibration_rejects_mismatched_unr_type(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    lookup_file, tenv8_dir = sample_unr_data
    with pytest.raises(ValueError, match="does not match the calibration grid_type"):
        run_calibration(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=lookup_file,
            unr_timeseries_dir=tenv8_dir,
            output_dir=tmp_path / "out",
            los_file=sample_static_los,
            calibration_reference_type="variable",
        )


def test_run_calibration_requires_staged_grid_type(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    lookup_file, tenv8_dir = sample_unr_data
    for path in tenv8_dir.glob("*_constant.tenv8"):
        path.unlink()
    with pytest.raises(FileNotFoundError, match="--grid-type constant"):
        run_calibration(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=lookup_file,
            unr_timeseries_dir=tenv8_dir,
            output_dir=tmp_path / "out",
            los_file=sample_static_los,
        )


@pytest.fixture
def captured_core(monkeypatch):
    """Replace Venti's calibrate_pair with a recorder returning a zero surface."""
    from venti.calibration import two_pass

    calls: list[dict] = []

    def _fake(disp, gnss_los, mask, ref_point, options, pixel_m, cycle_m, **kwargs):
        calls.append(
            {
                "disp": disp,
                "gnss_los": gnss_los,
                "mask": mask,
                "ref_point": ref_point,
                "options": options,
                "pixel_m": pixel_m,
                "cycle_m": cycle_m,
                **kwargs,
            }
        )
        zeros = np.zeros_like(disp)
        tropo = kwargs.get("tropo")
        set_corr = kwargs.get("set_correction")
        return two_pass.CalibrationResult(
            cal_gnss_surface=zeros,
            cal_reference_offset=zeros,
            cal_tropo=zeros if tropo is None else np.asarray(tropo, disp.dtype),
            cal_set=zeros if set_corr is None else np.asarray(set_corr, disp.dtype),
            cal_unwrap_shift=zeros,
            components_applied={
                "cal_gnss_surface": True,
                "cal_reference_offset": True,
                "cal_tropo": tropo is not None,
                "cal_set": set_corr is not None,
                "cal_unwrap_shift": False,
            },
            method=options.surface.method,
            reference_point=tuple(ref_point),
            reference_value=0.0,
        )

    monkeypatch.setattr(two_pass, "calibrate_pair", _fake)
    return calls


@pytest.mark.parametrize("apply_set", [True, False])
def test_run_calibration_venti_core_inputs(
    tmp_path: Path,
    sample_disp_product_with_corrections: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
    apply_set: bool,
):
    lookup_file, tenv8_dir = sample_unr_data
    options = gamma_parameters(
        {
            "calibration_options": {
                "apply_solid_earth_tide_correction": apply_set,
                "downsample_factor": 2,
                "downsample_weighted": True,
                "event_mask_buffer_pixels": 2,
            }
        }
    ).calibration_options

    out_path = run_calibration(
        disp_file=sample_disp_product_with_corrections,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        algorithm_parameters=AlgorithmParameters(calibration_options=options),
    )

    (call,) = captured_core
    with xr.open_dataset(
        sample_disp_product_with_corrections, group="corrections"
    ) as corr:
        set_layer = corr["solid_earth_tide"].values
    with xr.open_dataset(sample_disp_product_with_corrections) as disp:
        coherence = disp["temporal_coherence"].values

    # SET is handed to Venti as a component (removed before the fit, restored)
    if apply_set:
        np.testing.assert_allclose(call["set_correction"], set_layer, rtol=1e-6)
    else:
        assert call["set_correction"] is None
    assert call["tropo"] is None
    # Displacement in m, GNSS converted from mm to m, valid reference pixel
    assert call["disp"].dtype == np.float32
    assert np.nanmax(np.abs(call["gnss_los"])) < 0.01
    assert call["mask"][call["ref_point"]]
    assert round(call["options"].window_size_meters / call["pixel_m"]) == (
        options.window_size_pixels
    )
    assert call["options"].model_dump(
        include={"apply_solid_earth_tide_correction", "event_mask_buffer_pixels"}
    ) == {
        "apply_solid_earth_tide_correction": apply_set,
        "event_mask_buffer_pixels": 2,
    }
    assert call["options"].downsample_factor == 2
    assert call["options"].downsample_weighted
    np.testing.assert_allclose(call["coherence"], coherence, rtol=1e-6)
    assert call["options"].surface.method == "windowed_plane"  # gamma

    with xr.open_dataset(out_path) as ds:
        applied = ds["calibration"].attrs["corrections_applied"]
    assert applied == ("solid_earth_tide" if apply_set else "none")


def test_unwrap_error_correction_is_off_by_default():
    """Venti's mask-island unwrap correction is unvalidated: default off."""
    assert CalibrationOptions().unwrap_error_correction is False
    assert AlgorithmParameters().calibration_options.unwrap_error_correction is False
    assert CalibrationOptions().to_venti().unwrap_error_correction is False


@pytest.mark.parametrize("enabled", [True, False])
def test_run_calibration_passes_half_the_wavelength_as_the_cycle(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
    enabled: bool,
):
    """calibrate_pair takes the unwrapping cycle itself: lambda/2 of LOS
    displacement (two-way path), from the product's radar wavelength."""
    lookup_file, tenv8_dir = sample_unr_data
    run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        algorithm_parameters=AlgorithmParameters(
            calibration_options=CalibrationOptions(
                unwrap_error_correction=enabled,
                apply_solid_earth_tide_correction=False,
            )
        ),
    )

    (call,) = captured_core
    assert call["options"].unwrap_error_correction is enabled
    # The fixture product's /identification/radar_wavelength is 0.05546 m; one
    # unwrapping cycle of LOS displacement is lambda/2
    assert call["cycle_m"] == pytest.approx(0.05546 / 2, rel=1e-6)


def _displacement_with_nodata(tmp_path: Path) -> tuple[Path, np.ndarray]:
    """A float64 displacement file with NaN nodata, and its float64 values."""
    rng = np.random.default_rng(7)
    values = rng.normal(0.0, 0.05, (37, 23))  # not representable in float32
    values[3:9, 5:11] = np.nan
    values[20, :] = np.nan
    path = tmp_path / "displacement.nc"
    xr.Dataset({"displacement": (["y", "x"], values)}).to_netcdf(
        path, engine="h5netcdf"
    )
    return path, values


@pytest.mark.parametrize("block_rows", [1, 7, 512])
def test_load_displacement_matches_astype_float32(tmp_path: Path, block_rows: int):
    from cal_disp.workflow import _load_displacement

    path, values = _displacement_with_nodata(tmp_path)
    with xr.open_dataset(path, engine="h5netcdf") as ds:
        assert ds["displacement"].dtype == np.float64
        loaded = _load_displacement(ds, block_rows=block_rows)
        # Read in blocks: the float64 layer is not cached in the dataset
        assert not ds["displacement"].variable._in_memory

    assert loaded.dtype == np.float32
    # Bit-identical to the previous `displacement.values.astype(np.float32)`
    expected = values.astype(np.float32)
    assert loaded.tobytes() == expected.tobytes()
    # The nodata count from the float32 array equals the float64 one
    assert np.count_nonzero(np.isnan(loaded)) == int(np.isnan(values).sum()) == 59


def test_run_calibration_loads_displacement_once(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
):
    """The block-wise float32 load hands Venti the same array and mask as
    ``displacement.values.astype(np.float32)`` did."""
    import h5py

    with h5py.File(sample_disp_product, "a") as f:
        f["displacement"][10:15, 40:60] = np.nan  # 100 nodata pixels
    with xr.open_dataset(sample_disp_product) as disp:
        values = disp["displacement"].values
    assert values.dtype == np.float64

    lookup_file, tenv8_dir = sample_unr_data
    run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
    )

    (call,) = captured_core
    assert call["disp"].dtype == np.float32
    assert call["disp"].tobytes() == values.astype(np.float32).tobytes()
    np.testing.assert_array_equal(call["mask"], ~np.isnan(values))
    assert int((~call["mask"]).sum()) == 100


def test_run_calibration_without_set_layer_warns(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
    caplog,
):
    lookup_file, tenv8_dir = sample_unr_data
    run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
    )

    assert captured_core[0]["set_correction"] is None
    assert "calibrating without SET" in caplog.text


def test_find_reference_point():
    from cal_disp.workflow import _find_reference_point

    coherence = np.full((40, 50), 0.3, dtype=np.float32)
    coherence[5:15, 30:45] = 0.99  # best block
    coherence[30, 5] = 1.0  # higher, but masked out
    mask = np.ones(coherence.shape, dtype=bool)
    mask[30, 5] = False

    row, col = _find_reference_point(coherence, mask)

    assert 5 <= row < 15
    assert 30 <= col < 45
    # no coherence layer: every valid pixel ties, the rule still returns one
    r0, c0 = _find_reference_point(None, mask)
    assert mask[r0, c0]


@pytest.mark.parametrize("dtype", [bool, np.float32])
def test_load_layer_matches_values_and_caches_nothing(tmp_path: Path, dtype):
    """Masks and coherence are read in row blocks: same values as
    ``.values.astype``, and the dataset keeps no decoded copy."""
    from cal_disp.workflow import _load_layer

    rng = np.random.default_rng(1)
    values = rng.random((37, 23)).astype(np.float32)
    values[2, :] = np.nan
    if dtype is bool:
        values = (values > 0.5).astype(np.float32)
    path = tmp_path / "layer.nc"
    xr.Dataset({"layer": (["y", "x"], values)}).to_netcdf(path, engine="h5netcdf")
    with xr.open_dataset(path, engine="h5netcdf") as ds:
        loaded = _load_layer(ds, "layer", dtype, block_rows=5)
        assert not ds["layer"].variable._in_memory
        np.testing.assert_array_equal(loaded, ds["layer"].values.astype(dtype))
    assert loaded.dtype == np.dtype(dtype)


@pytest.mark.parametrize("seed", range(6))
def test_reference_rule_matches_opera_utils(tmp_path: Path, seed: int):
    """The ported rule picks the same pixel as opera_utils' find_reference_point."""
    rebase_reference = pytest.importorskip("opera_utils.disp.rebase_reference")
    from rasterio.transform import from_origin

    from cal_disp.workflow import _quality_reference_point

    rng = np.random.default_rng(seed)
    quality = rng.uniform(0.0, 1.0, (60, 80)).astype(np.float32)
    quality[rng.uniform(size=quality.shape) < 0.3] = 0.0  # masked pixels
    if seed == 5:
        quality = np.minimum(quality, 0.7)  # forces the lower thresholds
    path = tmp_path / "q.tif"
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=60,
        width=80,
        count=1,
        dtype="float32",
        transform=from_origin(0, 0, 30, 30),
    ) as dst:
        dst.write(quality, 1)
    expected = tuple(int(v) for v in rebase_reference.find_reference_point(path))
    assert _quality_reference_point(quality) == expected


def test_reference_rule_without_candidates_raises():
    from cal_disp.workflow import _quality_reference_point

    with pytest.raises(ValueError, match="No valid candidates"):
        _quality_reference_point(np.full((5, 5), 0.2, dtype=np.float32))


def test_run_calibration_clips_negative_std(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    monkeypatch,
):
    """RBF-interpolated GNSS sigmas can overshoot below 0: never written."""
    import venti.gnss

    def _overshooting_std(**kwargs):
        std = np.full(kwargs["los_east"].shape, 0.5, dtype=np.float32)  # mm
        std[: len(std) // 2] = -0.3
        return std

    monkeypatch.setattr(venti.gnss, "compute_gnss_los_std", _overshooting_std)
    lookup_file, tenv8_dir = sample_unr_data
    out_path = run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        algorithm_parameters=gamma_parameters(),  # the gamma route's RBF sigma
    )

    with xr.open_dataset(out_path) as ds:
        std = ds["calibration_std"].values[0]
    assert std.min() == 0  # negative sigmas clipped, not written
    np.testing.assert_array_equal(std[:80], 0)
    np.testing.assert_allclose(std[120:], 0.0005)


def test_gnss_field_on_fit_grid_matches_full_resolution(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
):
    """With downsample_factor > 1 the GNSS field is computed on the coarse
    fit grid and upsampled; it must match the full-resolution field."""
    import rasterio

    # Real LOS rasters vary smoothly; the fixture's is random per pixel, which
    # would change the LOS each station samples on the coarse grid.
    smooth_los = tmp_path / "smooth_los.tif"
    with rasterio.open(sample_static_los) as src:
        profile = src.profile
    xx = np.linspace(0, 1, 200, dtype=np.float32)[None, :].repeat(200, axis=0)
    east, north = -0.6 + 0.05 * xx, np.full_like(xx, -0.1)
    up = np.sqrt(1 - east**2 - north**2)
    with rasterio.open(smooth_los, "w", **profile) as dst:
        dst.write(np.stack([east, north, up]))

    lookup_file, tenv8_dir = sample_unr_data
    for factor in (1, 4):
        run_calibration(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=lookup_file,
            unr_timeseries_dir=tenv8_dir,
            output_dir=tmp_path / f"out{factor}",
            los_file=smooth_los,
            algorithm_parameters=AlgorithmParameters(
                calibration_options=CalibrationOptions(
                    grid_type="variable", downsample_factor=factor
                )
            ),
        )
    full, coarse = (call["gnss_los"] for call in captured_core)
    cache = np.load(
        tmp_path
        / "out4/scratch/gnss"
        / next(
            p.name
            for p in (tmp_path / "out4/scratch/gnss").glob("gnss_los_disp_*.npy")
            if "std" not in p.name
        )
    )

    assert cache.shape == (50, 50)  # computed on the fit grid
    assert coarse.shape == full.shape == (200, 200)
    spread = np.ptp(full)
    assert spread > 0
    # Between the outermost block centres (pixels 2..198) it is interpolated
    inner = (slice(2, 199), slice(2, 199))
    np.testing.assert_allclose(coarse[inner], full[inner], atol=0.005 * spread)
    # In the 2-pixel edge band the edge value is held
    np.testing.assert_allclose(coarse, full, atol=0.03 * spread)


@pytest.mark.parametrize("apply_tropo", [True, False])
def test_run_calibration_tropo_toggle(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
    monkeypatch,
    apply_tropo: bool,
):
    """apply_tropo_correction=False ignores the tropo files in the runconfig."""
    import rasterio

    from cal_disp.prep import tropo

    prepared = []

    def _fake_prepare(**_kwargs):
        # Differential delay of 1 cm everywhere (sec - ref)
        with rasterio.open(sample_static_los) as src:
            profile = src.profile | {"count": 1, "dtype": "float32"}
        paths = []
        for name, value in (("ref", 0.0), ("sec", 0.01)):
            path = tmp_path / f"{name}_tropo.tif"
            with rasterio.open(path, "w", **profile) as dst:
                dst.write(np.full((1, 200, 200), value, dtype=np.float32))
            paths.append(path)
        prepared.append(paths)
        return tuple(paths)

    monkeypatch.setattr(tropo, "prepare_troposphere_correction", _fake_prepare)
    lookup_file, tenv8_dir = sample_unr_data
    out_path = run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        dem_file=sample_static_dem,
        reference_tropo_files=[tmp_path / "ref.nc"],
        secondary_tropo_files=[tmp_path / "sec.nc"],
        algorithm_parameters=AlgorithmParameters(
            calibration_options=CalibrationOptions(apply_tropo_correction=apply_tropo)
        ),
    )

    tropo_component = captured_core[0]["tropo"]
    with xr.open_dataset(out_path) as ds:
        applied = ds["calibration"].attrs["corrections_applied"]
    with xr.open_dataset(out_path, group="metadata") as meta:
        atmospheric = str(meta["ceos_atmospheric_phase_correction"].values)
    if apply_tropo:
        assert len(prepared) == 1
        np.testing.assert_allclose(tropo_component, 0.01)
        assert applied == "troposphere"
        assert atmospheric == "tropospheric"
    else:
        assert prepared == []
        assert tropo_component is None
        assert applied == "none"
        assert atmospheric == "none"


def _write_mask_like(disp_file: Path, path: Path, values: np.ndarray, **overrides):
    """Write a uint8 mask raster on the grid of `disp_file` (or a shifted one)."""
    from rasterio.transform import from_origin

    with xr.open_dataset(disp_file, engine="h5netcdf") as ds:
        x, y = ds.x.values, ds.y.values
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    profile = {
        "driver": "GTiff",
        "height": values.shape[0],
        "width": values.shape[1],
        "count": 1,
        "dtype": values.dtype,
        "transform": from_origin(
            float(x[0]) - dx / 2 + overrides.pop("shift_x", 0.0),
            float(y[0]) - dy / 2,
            dx,
            -dy,
        ),
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(values, 1)
    return path


def _calibrate(tmp_path, disp, los, unr, **kwargs):
    lookup_file, tenv8_dir = unr
    return run_calibration(
        disp_file=disp,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=los,
        algorithm_parameters=AlgorithmParameters(
            calibration_options=CalibrationOptions(
                apply_solid_earth_tide_correction=False
            )
        ),
        **kwargs,
    )


def test_mask_file_pixels_are_invalid_for_the_fit(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
):
    """The runconfig mask_file is ANDed into the valid mask Venti fits on."""
    _calibrate(tmp_path / "a", sample_disp_product, sample_static_los, sample_unr_data)
    (plain,) = captured_core
    captured_core.clear()
    values = np.ones(plain["mask"].shape, dtype=np.uint8)
    values[50:80, 60:120] = 0
    mask_file = _write_mask_like(sample_disp_product, tmp_path / "mask.tif", values)
    _calibrate(
        tmp_path / "b",
        sample_disp_product,
        sample_static_los,
        sample_unr_data,
        mask_file=mask_file,
    )
    (masked,) = captured_core
    assert not masked["mask"][50:80, 60:120].any()
    outside = np.ones(values.shape, bool)
    outside[50:80, 60:120] = False
    np.testing.assert_array_equal(masked["mask"][outside], plain["mask"][outside])


def test_frame_fit_mask_is_read_from_the_package_and_passed_to_the_fit(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
    monkeypatch,
):
    """surface.fit_mask (frame table) names a packaged mask: it reaches
    calibrate_pair as fit_mask, and the valid mask itself is unchanged."""
    import venti.frames

    _calibrate(tmp_path / "a", sample_disp_product, sample_static_los, sample_unr_data)
    (plain,) = captured_core
    captured_core.clear()
    assert plain["fit_mask"] is None
    values = np.ones(plain["mask"].shape, dtype=np.uint8)
    values[50:80, 60:120] = 0
    masks = tmp_path / "fit_masks"
    masks.mkdir()
    _write_mask_like(sample_disp_product, masks / "F00000_stable.tif", values)
    monkeypatch.setattr(venti.frames, "FIT_MASK_DIR", masks)
    lookup_file, tenv8_dir = sample_unr_data
    run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "b",
        los_file=sample_static_los,
        algorithm_parameters=AlgorithmParameters(
            calibration_options=CalibrationOptions(
                apply_solid_earth_tide_correction=False,
                surface={"fit_mask": "F00000_stable.tif"},
            )
        ),
    )
    (fitted,) = captured_core
    np.testing.assert_array_equal(fitted["fit_mask"], values == 1)
    np.testing.assert_array_equal(fitted["mask"], plain["mask"])


@pytest.mark.usefixtures("captured_core")
def test_unknown_fit_mask_name_fails_before_the_fit(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    lookup_file, tenv8_dir = sample_unr_data
    with pytest.raises(FileNotFoundError, match="no packaged fit mask"):
        run_calibration(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=lookup_file,
            unr_timeseries_dir=tenv8_dir,
            output_dir=tmp_path / "out",
            los_file=sample_static_los,
            algorithm_parameters=AlgorithmParameters(
                calibration_options=CalibrationOptions(
                    apply_solid_earth_tide_correction=False,
                    surface={"fit_mask": "F00000_missing.tif"},
                )
            ),
        )


@pytest.mark.usefixtures("captured_core")  # Venti stubbed: the error comes first
@pytest.mark.parametrize(
    ("problem", "match"),
    [
        ("dtype", "dtype must be uint8"),
        ("shape", "does not match the DISP grid"),
        ("grid", "is not the DISP grid"),
        ("values", "values must be 0 or 1"),
    ],
)
def test_mask_file_must_be_a_0_1_uint8_raster_on_the_disp_grid(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    problem: str,
    match: str,
):
    with xr.open_dataset(sample_disp_product, engine="h5netcdf") as ds:
        shape = (ds.sizes["y"], ds.sizes["x"])
    values = np.ones(shape, dtype=np.uint8)
    extra = {}
    if problem == "dtype":
        values = values.astype(np.float32)
    elif problem == "shape":
        values = values[:-3]
    elif problem == "grid":
        extra["shift_x"] = 500.0
    elif problem == "values":
        values[0, 0] = 255
    mask_file = _write_mask_like(
        sample_disp_product, tmp_path / "m.tif", values, **extra
    )
    with pytest.raises(ValueError, match=match):
        _calibrate(
            tmp_path,
            sample_disp_product,
            sample_static_los,
            sample_unr_data,
            mask_file=mask_file,
        )


# --- Venti calibrate_pair wiring (Venti plan T37) ----------------------------


def _fake_tropo(monkeypatch, tmp_path, los_file, sec_field):
    """prepare_troposphere_correction writing ref = 0 and sec = `sec_field`."""
    from cal_disp.prep import tropo

    prepared = []

    def _fake_prepare(**_kwargs):
        with rasterio.open(los_file) as src:
            profile = src.profile | {"count": 1, "dtype": "float32"}
        paths = []
        for name, field in (("ref", np.zeros((200, 200))), ("sec", sec_field)):
            path = tmp_path / f"{name}_tropo.tif"
            with rasterio.open(path, "w", **profile) as dst:
                dst.write(field[np.newaxis].astype(np.float32))
            paths.append(path)
        prepared.append(paths)
        return tuple(paths)

    monkeypatch.setattr(tropo, "prepare_troposphere_correction", _fake_prepare)
    return prepared


@pytest.mark.parametrize("mode", ["off", "full", "stratified"])
def test_tropo_modes_reach_calibrate_pair(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_static_dem: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
    monkeypatch,
    mode: str,
):
    """tropo.mode: off ignores the files; full passes sec - ref; stratified
    passes Venti's height-dependent part of it (Venti plan T34 through T37)."""
    from venti.calibration.tropo import apply_tropo

    with rasterio.open(sample_static_dem) as src:
        dem = src.read(1).astype(np.float64)
    yy, xx = np.mgrid[0:200, 0:200]
    # 4 mm/km stratification plus a lateral pattern the model cannot explain
    sec = -0.004 * dem / 1000.0 + 0.003 * np.sin(xx / 7.0) * np.cos(yy / 9.0)
    prepared = _fake_tropo(monkeypatch, tmp_path, sample_static_los, sec)
    lookup_file, tenv8_dir = sample_unr_data
    out_path = run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        dem_file=sample_static_dem,
        reference_tropo_files=[tmp_path / "ref.nc"],
        secondary_tropo_files=[tmp_path / "sec.nc"],
        algorithm_parameters=AlgorithmParameters.create_default().with_overrides(
            {"calibration_options": {"tropo": {"mode": mode}}}
        ),
    )
    (call,) = captured_core
    with xr.open_dataset(out_path) as ds:
        applied = ds["calibration"].attrs["corrections_applied"]
    if mode == "off":
        assert prepared == []
        assert call["tropo"] is None
        assert "troposphere" not in applied
        return
    assert len(prepared) == 1
    assert "troposphere" in applied
    full = sec.astype(np.float32)
    if mode == "full":
        np.testing.assert_allclose(call["tropo"], full, atol=1e-7)
    else:
        expected = apply_tropo(
            full,
            "stratified",
            dem=dem,
            valid=call["mask"],
            x=np.asarray(xr.open_dataset(sample_disp_product).x),
            y=np.asarray(xr.open_dataset(sample_disp_product).y),
        )
        np.testing.assert_allclose(call["tropo"], expected, atol=1e-6)
        # the lateral pattern is not passed on
        assert np.nanstd(call["tropo"] - (-0.004 * dem / 1000.0)) < 0.0005


def test_loclin_end_to_end_writes_sigma_cal(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    """surface.method = loclin through the real Venti: the product carries the
    loclin calibration and sigma_CAL as calibration_std."""
    lookup_file, tenv8_dir = sample_unr_data
    params = AlgorithmParameters.create_default().with_overrides(
        {
            "calibration_options": {
                "apply_solid_earth_tide_correction": False,
                "downsample_factor": 2,
                "surface": {"method": "loclin", "cutoff_wavelength_meters": 5000.0},
                "uncertainty": {"k_grid": 3.9},
            }
        }
    )
    out = run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        algorithm_parameters=params,
    )
    with xr.open_dataset(out) as ds:
        cal = ds["calibration"].values
        std = ds["calibration_std"].values
        source = ds["calibration_std"].attrs["uncertainty_source"]
    assert np.isfinite(cal).mean() > 0.9
    assert "sigma_CAL" in source and "k = 3.9" in source
    assert np.nanmin(std) > 0


def test_gnss_provenance_written_on_v05_route_only(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    """T37.3: the loclin (v0.5) route records the GNSS realisation in
    /metadata/gnss_provenance, readable by VLM (venti.gnss.vlm); the gamma
    route writes none, so the gamma golden keeps its structure."""
    from venti.gnss.vlm import PROVENANCE_KEYS, read_cal_provenance

    lookup_file, tenv8_dir = sample_unr_data

    def run(overrides: dict, name: str) -> Path:
        params = AlgorithmParameters.create_default().with_overrides(
            {
                "calibration_options": {
                    "apply_solid_earth_tide_correction": False,
                    "downsample_factor": 2,
                    **overrides,
                }
            }
        )
        return run_calibration(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=lookup_file,
            unr_timeseries_dir=tenv8_dir,
            output_dir=tmp_path / name,
            los_file=sample_static_los,
            algorithm_parameters=params,
        )

    v05 = run(
        {"surface": {"method": "loclin", "cutoff_wavelength_meters": 5000.0}}, "v05"
    )
    record = read_cal_provenance(v05)
    assert record is not None
    assert set(PROVENANCE_KEYS) <= set(record)
    assert record["grid_type"] == "constant" and len(record["digest"]) == 64
    assert (
        read_cal_provenance(run({"surface": {"method": "windowed_plane"}}, "gamma"))
        is None
    )


def test_algorithm_options_reach_venti():
    """The schema-v2 groups of algorithm_parameters arrive in Venti unchanged."""
    params = AlgorithmParameters.create_default().with_overrides(
        {
            "calibration_options": {
                "downsample_factor": 3,
                "surface": {"method": "loclin", "two_pass": True},
                "weights": {"coherence_power": 8.0, "robust": True},
                "tropo": {"mode": "auto"},
                "uncertainty": {"k_grid": 3.9},
            }
        }
    )
    venti = params.calibration_options.to_venti()
    assert venti.downsample_factor == 3
    assert (venti.surface.method, venti.surface.two_pass) == ("loclin", True)
    assert (venti.weights.coherence_power, venti.weights.robust) == (8.0, True)
    assert venti.tropo.mode == "auto"
    assert venti.uncertainty.k_grid == 3.9
    default = AlgorithmParameters.create_default().calibration_options.to_venti()
    assert default.surface.method == "loclin"  # v0.5 by default since T49
    assert default.tropo.mode == "legacy"


@pytest.mark.parametrize(
    ("method", "enabled", "expect_hook"),
    [("loclin", True, True), ("loclin", False, False), ("windowed_plane", True, False)],
)
def test_unwrap_hook_only_on_the_loclin_route(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
    method: str,
    enabled: bool,
    expect_hook: bool,
):
    """Regression: with loclin, unwrap_error_correction=True used to do nothing
    because no hook reached calibrate_pair (the gamma route has its own)."""
    lookup_file, tenv8_dir = sample_unr_data
    params = AlgorithmParameters.create_default().with_overrides(
        {
            "calibration_options": {
                "apply_solid_earth_tide_correction": False,
                "unwrap_error_correction": enabled,
                "downsample_factor": 2,
                "surface": {"method": method, "cutoff_wavelength_meters": 5000.0},
            }
        }
    )
    run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        algorithm_parameters=params,
    )
    (call,) = captured_core
    assert callable(call["unwrap_hook"]) is expect_hook


def test_loclin_unwrap_end_to_end_writes_decisions(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    """loclin + unwrap through the real Venti: the hook runs and its per-region
    decisions are written next to the product."""
    lookup_file, tenv8_dir = sample_unr_data
    params = AlgorithmParameters.create_default().with_overrides(
        {
            "calibration_options": {
                "apply_solid_earth_tide_correction": False,
                "unwrap_error_correction": True,
                "downsample_factor": 2,
                "surface": {
                    "method": "loclin",
                    "cutoff_wavelength_meters": 5000.0,
                    "two_pass": True,
                },
            }
        }
    )
    out = run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        algorithm_parameters=params,
    )
    assert (tmp_path / "out" / "unwrap_decisions.csv").exists()
    with xr.open_dataset(out) as ds:
        assert np.isfinite(ds["calibration"].values).mean() > 0.9


@pytest.mark.parametrize("s0", [0.0, 0.003])
def test_calibration_std_states_what_sigma_cal_covers(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    s0: float,
):
    """sigma_CAL is the calibration-surface uncertainty (owner, 2026-10-09);
    the optional non-secular term is named with its value when set, and the
    users are told to add it when not."""
    lookup_file, tenv8_dir = sample_unr_data
    params = AlgorithmParameters.create_default().with_overrides(
        {
            "calibration_options": {
                "apply_solid_earth_tide_correction": False,
                "downsample_factor": 2,
                "surface": {"method": "loclin", "cutoff_wavelength_meters": 5000.0},
                "uncertainty": {"k_grid": 3.9, "sigma_nonsecular_meters": s0},
            }
        }
    )
    out = run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        algorithm_parameters=params,
    )
    with xr.open_dataset(out) as ds:
        source = ds["calibration_std"].attrs["uncertainty_source"]
    assert "uncertainty of the calibration surface" in source
    if s0:
        assert "sigma_nonsecular = 3 mm" in source
    else:
        assert "non-secular station motion" in source


@pytest.mark.parametrize(
    ("method", "expect_enu"), [("loclin", True), ("windowed_plane", False)]
)
def test_gnss_buffer_reaches_the_sampler_on_the_loclin_route(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    monkeypatch,
    method: str,
    expect_enu: bool,
):
    """v0.5 samples the grid with sample_gnss_enu and the gnss.buffer_meters
    option (R-G2); the gamma route keeps its fit path (golden unchanged).

    Both routes also sample the coarse `/auxiliary` grid; only calls on the
    fit grid count here.
    """
    import venti.gnss.sampling as sampling

    from cal_disp.workflow import AUXILIARY_STEP

    with xr.open_dataset(sample_disp_product) as ds:
        n_aux = len(ds.x.values[::AUXILIARY_STEP])
    seen: list = []
    aux: list = []
    real = sampling.sample_gnss_enu

    def spy(cfg, grid, exclude=None, defo_db_version=None):
        (aux if len(grid[0]) == n_aux else seen).append(cfg)
        return real(cfg, grid, exclude, defo_db_version)

    monkeypatch.setattr(sampling, "sample_gnss_enu", spy)
    lookup_file, tenv8_dir = sample_unr_data
    params = AlgorithmParameters.create_default().with_overrides(
        {
            "calibration_options": {
                "apply_solid_earth_tide_correction": False,
                "downsample_factor": 2,
                "surface": {"method": method, "cutoff_wavelength_meters": 5000.0},
                "gnss": {"buffer_meters": 20000.0},
            }
        }
    )
    out = run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        algorithm_parameters=params,
    )
    assert bool(seen) is expect_enu
    assert len(aux) == 1 and aux[0].buffer_meters == 20000.0
    if expect_enu:
        assert seen[0].buffer_meters == 20000.0
        with xr.open_dataset(out) as ds:
            assert np.isfinite(ds["calibration"].values).mean() > 0.9


def test_loclin_remove_restore_reads_the_venti_databases(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
    captured_core,
):
    """v0.5: areas always excluded, events only for pairs that span them,
    from the Venti GeoJSON (id/name/version, t0) on the product grid."""
    import json

    from pyproj import Transformer

    to_ll = Transformer.from_crs("EPSG:32611", "EPSG:4326", always_xy=True)

    def box(x0, x1, y0, y1):
        ring = [
            to_ll.transform(x, y)
            for x, y in [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
        ]
        return {"type": "Polygon", "coordinates": [ring]}

    # product grid: x 405000-424900, y 3778000 down to 3758100 (100 m)
    areas = {
        "type": "FeatureCollection",
        "version": "t-1",
        "features": [
            {
                "type": "Feature",
                "properties": {"id": "west", "name": "west half"},
                "geometry": box(404000, 414950, 3757000, 3779000),
            }
        ],
    }
    events = {
        "type": "FeatureCollection",
        "version": "t-1",
        "features": [
            {
                "type": "Feature",
                "properties": {"id": "in_pair", "t0": "2022-03-01"},
                "geometry": box(415050, 420000, 3757000, 3779000),
            },
            {
                "type": "Feature",
                "properties": {"id": "later", "t0": "2023-03-01"},
                "geometry": box(420050, 426000, 3757000, 3779000),
            },
        ],
    }
    (tmp_path / "areas.geojson").write_text(json.dumps(areas))
    (tmp_path / "events.geojson").write_text(json.dumps(events))
    lookup_file, tenv8_dir = sample_unr_data
    params = AlgorithmParameters.create_default().with_overrides(
        {
            "calibration_options": {
                "apply_solid_earth_tide_correction": False,
                "downsample_factor": 2,
                "surface": {"method": "loclin", "cutoff_wavelength_meters": 5000.0},
                # the fixture has four grid nodes: excluding those in the area
                # leaves too few to re-interpolate from (the v0.5 default)
                "gnss": {"exclude_defo_nodes": False, "reinterpolate_excluded": False},
            }
        }
    )
    run_calibration(
        disp_file=sample_disp_product,
        unr_grid_latlon_file=lookup_file,
        unr_timeseries_dir=tenv8_dir,
        output_dir=tmp_path / "out",
        los_file=sample_static_los,
        algorithm_parameters=params,
        defo_area_db_json=tmp_path / "areas.geojson",
        event_db_json=tmp_path / "events.geojson",
    )
    (call,) = captured_core
    exclude = call["exclude_mask"]
    assert exclude[:, :95].all()  # the area, always
    assert exclude[:, 102:148].all()  # the event inside the pair
    assert not exclude[:, 153:].any()  # the later event does not apply


def test_v05_metadata_records_schema_frame_table_and_components(
    tmp_path: Path,
    sample_disp_product: Path,
    sample_static_los: Path,
    sample_unr_data: tuple[Path, Path],
):
    """T37.3 / T27.4: the v0.5 product records the algorithm schema, the
    frame-parameter table version and the applied components; gamma does not."""
    import json

    import h5py

    lookup_file, tenv8_dir = sample_unr_data

    def run(params, name, **kw):
        return run_calibration(
            disp_file=sample_disp_product,
            unr_grid_latlon_file=lookup_file,
            unr_timeseries_dir=tenv8_dir,
            output_dir=tmp_path / name,
            los_file=sample_static_los,
            algorithm_parameters=params,
            **kw,
        )

    v05 = run(
        AlgorithmParameters.create_default().with_overrides(
            {
                "calibration_options": {
                    "apply_solid_earth_tide_correction": False,
                    "downsample_factor": 2,
                    "surface": {"cutoff_wavelength_meters": 5000.0},
                }
            }
        ),
        "v05",
        frame_parameters_version="0.6-draft",
    )
    with h5py.File(v05) as f:
        meta = f["metadata"]
        assert int(meta["algorithm_schema_version"][()]) == 2
        assert meta["frame_parameters_version"][()].decode() == "0.6-draft"
        comps = json.loads(meta["components_applied"][()])
    assert comps["cal_gnss_surface"] is True
    assert comps["cal_set"] is False  # SET off in this run
    gamma = run(gamma_parameters(), "gamma", frame_parameters_version="0.6-draft")
    with h5py.File(gamma) as f:
        for name in (
            "algorithm_schema_version",
            "frame_parameters_version",
            "components_applied",
        ):
            assert name not in f["metadata"]  # gamma keeps its 0.3 structure
