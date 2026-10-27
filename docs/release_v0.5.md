# Releasing v0.5 (Venti plan T49–T51)

The steps that turn the validated v0.5 algorithm into the default. Each one
changes what a run produces, so they happen together, after the VnV sign-off
(T50.3), in this order.

## 0. Before you start
- [x] VnV report signed off (owner 2026-10-15; OPERA CalVal pending) (disp2vlm_validation `docs/vnv/v0.5-calval-report.md`).
- [x] Decision on `calibration_std` NaN where the fit coverage is low (2026-10-14: extrapolated σ, `uncertainty.extrapolate_sigma`)
      (CHANGELOG, "Open before the release").
- [x] Frame-key fix merged (#16), so the frame table reaches runs.

## 1. Flip the defaults (T49.1)
- [x] `cp configs/algorithm_parameters_v05.yaml configs/algorithm_parameters.yaml`
      (keep the field documentation; only the values change).
- [x] Model defaults in `src/cal_disp/config/_algorithm.py` follow the file.
- [x] `scripts/build_golden_output.sh` writes the v0.5 parameters and adds the
      static ancillaries production uses:
      `algorithm_parameters_overrides_json` (Venti
      `load_frame_table().write(path, materialize=True)`), `defo_area_db_json`
      and `event_db_json` (Venti `data/`).
- [x] `tests/test_release_configs.py`: the gamma test reads its own copy of the
      0.3 parameters (the build script now writes v0.5).

## 2. Image (T49.2)
- [x] `docker build --network=host -t cal-disp:0.5.0-rc -f docker/Dockerfile .`
      (docker-build skill); the env budget check (`budget.yaml`) passes.

## 3. Golden (T49.3)
- [x] Regenerate inside the image: `scripts/build_golden_output.sh`.
- [x] `scripts/run_validation.sh --golden-dir test_golden` passes in the image.
- [x] New `MANIFEST.sha256`; archive the 0.3 golden under `test_golden/archive/`.
- [x] Gamma reproduction: a run with `configs/algorithm_parameters_gamma.yaml`
      matches the archived 0.3 golden.

## 4. Changelog (T49.4)
- [x] Re-measure the golden-pair table in the image and update `CHANGELOG.md`.
- [x] Date the `[0.5.0]` entry.

## 5. Release (T51)
- [x] Tag `v0.5.0` on the merge commit (`438b84e`, 2026-10-15); image `cal-disp:0.5.0` built from the tag and validated against the golden.
- [ ] Push the image to the registry (owner).

## 0.5.1 (2026-10-16)
- [x] Metadata fields (T37.3, T27.4; #24); golden rebuilt in `0.5.1-rc` (metadata only).
- [x] Tag `v0.5.1` on `54e6382`; image `cal-disp:0.5.1` built from the tag and validated against the golden.
- [ ] Push the image (owner, with 0.5.0's).

## 0.5.2 (2026-10-20)
- [x] Venti `ca58158` (frame table 0.7-draft: F18904 mirrors LA; venti-dev #44); golden rebuilt in `0.5.2-rc` (metadata only: `frame_parameters_version`).
- [x] Tag `v0.5.2` on `1d64f2f` (owner, 2026-10-20); image `cal-disp:0.5.2` built from the tag and validated against the golden (`test_golden/validation_0.5.2_release_image.log`).
- [ ] Push the image (owner). The local 0.5.0/0.5.1 images were pruned on 2026-10-20; rebuild from the tags if they are to be pushed.

## 0.5.3 (2026-10-20)
- [x] Venti `50689b6` (venti-dev #45: stratified tropo evaluated wherever DEM and delay exist); image `0.5.3-rc` validates against the 0.5.2 golden (no rebuild).
- [x] Tag `v0.5.3` on `3f2c241` (owner, 2026-10-20); image `cal-disp:0.5.3` built from the tag and validated against the golden (`test_golden/validation_0.5.3_release_image.log`).
- [ ] Push the image (owner).

## 0.5.4 (2026-10-22)
- [x] Venti `1fa4436` (frame table 0.8-draft: Hawaii stratified tropo, venti-dev #48; decision 0004); golden rebuilt in `0.5.4-rc` (metadata only: `frame_parameters_version`).
- [x] Tag `v0.5.4` on `c3e96bb` (owner, 2026-10-22); image `cal-disp:0.5.4` built from the tag and validated against the golden (`test_golden/validation_0.5.4_release_image.log`).
- [ ] Push the image (owner).

## 0.5.5 (2026-10-22)
- [x] `/plate_motion` group (ITRF2020-PMM, per-pair LOS) and GNSS E/N/U displacement in `/auxiliary` (#33); peak-memory fix (#34). `calibration` and `calibration_std` unchanged.
- [x] Golden rebuilt in `0.5.5-rc` (new `/plate_motion`, filled `/auxiliary`); reproduced in the image, peak 5.70 GB, budget OK.
- [x] Tag `v0.5.5` on `c1f3579` (owner, 2026-10-22; first placed on `8e5e797`, before #34, and moved); image `cal-disp:0.5.5` built from the tag and validated against the golden, peak 5.77 GB, budget OK (`test_golden/validation_0.5.5_release_image.log`).
- [ ] Push the image (owner).
- [ ] cal-vlm-batch: rebuild the cal-disp worker image on `cal-disp:0.5.5` so Batch products carry the new layers.

## 0.5.6 (2026-10-25)
- [x] Runconfig `fit_mask_file` (#37) on Venti `f2db0cb` (venti-dev #49, #50); no value change.
- [x] Golden unchanged: the 0.5.5 golden validates in `0.5.6-rc` on every group, peak 5.72 GB, budget OK.
- [x] Tag `v0.5.6` on `8d96246` (owner, 2026-10-25); image `cal-disp:0.5.6` built from the tag and validated against the golden, peak 5.83 GB, budget OK (`test_golden/validation_0.5.6_release_image.log`).
- [ ] Push the image (owner).
- [ ] cal-vlm-batch: worker on `cal-disp:0.5.6`; Gulf/East rerun with `--fit-mask-dir campaigns/gulf_east/fit_masks`.
