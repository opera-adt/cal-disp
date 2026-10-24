# Changelog

Product-affecting entries are marked **[product]**. Numbers on the golden pair
(F08882, 2022-01-11 → 2022-07-22) compare against the 0.3 golden.

## [0.5.5] — 2026-10-22 (tag `v0.5.5` on `c1f3579`)

### Added
- **[product]** `/plate_motion` group: LOS displacement of the frame's rigid
  plate rotation (ITRF2020-PMM, plate from the frame table, default NA) over
  the pair, on the DISP grid (`geepers.euler.plate_velocity_enu`, rounded to
  1 µm). `DISP − calibration` stays in IGS20; `DISP − calibration −
  plate_motion` is plate-fixed. Group attributes: `plate`,
  `plate_motion_model`, `pair_span_years`. +6 MB on the golden pair.
- `cal-disp validate` compares the `plate_motion` group (with `auxiliary`
  under `--group auxiliary`; absent from both products is a warning).

### Changed
- **[product]** `/auxiliary` `east_west`, `north_south`, `up_down` are no
  longer a zero placeholder: they are the GNSS E/N/U displacement over the
  pair from the UNR constant grid the calibration used (sampled with Venti and
  the same buffer/exclusion options), with new `*_std` layers; the group
  records `pair_span_years` (velocity = displacement / span). Variable grid:
  zeros as before, flagged in the group description.

### Fixed
- `plate_motion` is computed after the fit from the LOS file in row windows
  (#34): computed before the fit, it held one float32 frame through it and the
  golden run peaked at 6.02 GB, over the 6.0 GB budget. Now 5.70 GB in
  `cal-disp:0.5.5-rc` (0.5.4: 5.72–5.75 GB).

### Golden
Rebuilt in `cal-disp:0.5.5-rc` (`docs/golden/MANIFEST_0.5.5.sha256`).
`calibration`, `calibration_std` and the identification/metadata groups pass
against the 0.5.4 golden (`cal-disp validate --group main`); `/auxiliary`
changes and `/plate_motion` is new (+6 MB). Reproduced in the image at 1e-6,
peak 5.70 GB, dependency budget OK. The 0.5.4 golden is archived in
`test_golden/archive/golden_v0.5.4_20261022/`.

## [0.5.4] — 2026-10-22 (tag `v0.5.4` on `c3e96bb`)

### Changed
- **[product]** Venti `1fa4436`: frame table 0.8-draft gives Hawaii (F23211)
  `tropo.mode = stratified` (venti-dev #48). Stratified vs off on v0.5.3:
  DD sill −27.6% [−40.5, −19.4], velocity RMSE 2.40 → 2.30 mm/yr; the
  |median bias| rise on 18 stations is not significant (disp2vlm_validation
  decision 0004, accepted 2026-10-21). Values change only for F23211, which
  needs staged TROPO-ZENITH files.
- Venti also adds `surface.region_model` (`offset` default, unchanged values;
  `plane` opt-in for large islands, venti-dev #47).

### Golden
Rebuilt in `cal-disp:0.5.4-rc` (`docs/golden/MANIFEST_0.5.4.sha256`). Every data
layer is bit-identical to 0.5.2; the one product difference is
`/metadata/frame_parameters_version` (`0.7-draft` → `0.8-draft`), which made
the 0.5.2 golden fail (1 difference). `algorithm_parameters_yaml` gains the two
new options at their defaults (compared as parameters, equal). The 0.5.2 golden
is archived in `test_golden/archive/golden_v0.5.2_20261020/`.

## [0.5.3] — 2026-10-20 (tag `v0.5.3` on `3f2c241`)

### Fixed
- **[product]** Venti `50689b6` (venti-dev #45): the stratified tropo model is
  evaluated wherever the DEM and the delay exist; the fit still uses only the
  pair's valid pixels. Before, `cal_tropo` (and so `calibration`) was NaN
  wherever that epoch's DISP was masked, 16–24% of each epoch on the frames
  with `tropo.mode = stratified` (F16940, F18904, F35991), ~90% of it pixels
  with finite DISP and `recommended_mask = 0`. Pixels that were valid are
  bit-identical. LA VLM: 172 gated stations, NMAD 1.29 mm/yr (161, 1.32);
  Puerto Rico DISP gate PASS, velocity RMSE 3.18 → 3.08 mm/yr.

### Golden
Unchanged: the golden frame (F08882) runs with tropo off, and the 0.5.3 image
validates against the 0.5.2 golden (`test_golden/validation_0.5.3rc1.log`).

## [0.5.2] — 2026-10-20 (tag `v0.5.2` on `1d64f2f`)

### Changed
- **[product]** Venti `ca58158`: frame table 0.7-draft adds F18904, the LA
  descending VLM partner, with LA's settings (`surface.regions = separate`,
  `tropo.mode = stratified`; venti-dev #44). Unlisted, F18904's Catalina island
  offset was wrong (CAT2 LOS −54.2 vs GNSS −29.7 mm/yr). With the entry: F18904
  DISP gate PASS, DD sill −91.4% vs raw (−89.0% unlisted); LA VLM gated NMAD
  1.32 mm/yr (1.46). Values change only for F18904; F18904 needs staged
  TROPO-ZENITH files for the stratified correction.

### Golden
Rebuilt in `cal-disp:0.5.2-rc` (`docs/golden/MANIFEST_0.5.2.sha256`). Every data
layer is bit-identical to 0.5.1; the only product difference is
`/metadata/frame_parameters_version` (`0.6-draft` → `0.7-draft`), which made
the 0.5.1 golden fail (1 difference). The 0.5.1 golden is archived in
`test_golden/archive/golden_v0.5.1_20261016/`.

## [0.5.1] — 2026-10-16 (tag `v0.5.1` on `54e6382`)

### Added
- **[product]** `/metadata` on the v0.5 route records `algorithm_schema_version`
  (Venti schema, 2), `frame_parameters_version` (top-level version of the
  runconfig's frame table, e.g. `0.6-draft`) and `components_applied` (which
  `cal_*` components are summed into `calibration`, JSON). Venti plan T37.3 and
  T27.4. The gamma route is unchanged.

### Golden
Rebuilt in `cal-disp:0.5.1-rc` (`docs/golden/MANIFEST_0.5.1.sha256`). Every data
layer is bit-identical to 0.5.0; only the three metadata variables are new. The
0.5.0 golden is archived in `test_golden/archive/golden_v0.5.0_20261015/`.

### Fixed since 0.5.0
- Venti `35ec67e`: falls back to dropping excluded GNSS nodes when too few
  remain to re-interpolate (#22).
- `cal-disp validate` compares `gnss_provenance` field by field; the field
  hashes only warn, so validation also passes outside the image (#23).

## [0.5.0] — 2026-10-15 (Venti plan T49–T51; tag `v0.5.0` on `438b84e`)

v0.5 replaces the gamma windowed-plane surface with Venti's local-linear
surface and adds remove-restore, the 50 km GNSS buffer, per-frame settings and
σ_CAL. Validated on 10 benchmark frames (T46: all PASS under decisions 0002 and
0003; disp2vlm_validation `docs/vnv/v0.5-calval-report.md`), confirmed on all
10 frames with the release candidate through the frame table (2026-10-15).
VnV signed off by the owner 2026-10-15.

**The defaults are now v0.5 (T49.1).** `configs/algorithm_parameters.yaml` and
the model defaults are the v0.5 set. A parameter file that leaves out the
schema-v2 groups (`surface`, `weights`, `gnss`, `uncertainty`), as every 0.3
file does, now runs v0.5 for them: to reproduce gamma use
`configs/algorithm_parameters_gamma.yaml`, which spells them out.

### Added
- `configs/algorithm_parameters_v05.yaml` (the T46 set) and
  `configs/algorithm_parameters_gamma.yaml` (decision 0001, reproduces the 0.3
  golden); `tests/test_release_configs.py` keeps both in step.
- **[product]** v0.5 route through Venti `calibrate_pair` with schema-v2 options:
  `surface` (loclin, 50 km half-response cutoff, robust two-pass tie), `weights`
  (coherence⁸ × robust), `tropo` (per-frame mode), `uncertainty` (σ_CAL),
  `gnss`, `unwrap` (#6).
- **[product]** GNSS field from the UNR grid with nodes up to 50 km outside the
  frame (R-G2, R-G5) (#10).
- **[product]** Remove-restore from the Venti area and event databases
  (`defo_area_db_json`, `event_db_json`): deforming areas and event footprints
  are left out of the fit and restored; GNSS nodes inside them are excluded and
  re-interpolated (#11).
- **[product]** Per-frame stable-pixel fit mask (`surface.fit_mask`, packaged
  with Venti; F08882 only) (#15).
- **[product]** `/metadata/gnss_provenance` record on the v0.5 route (#13).
- Runconfig `mask_file` is applied (#4).
- Unwrap-error correction runs on the loclin route when enabled; off by
  default (owner 2026-10-09) (#7).
- Optional non-secular σ term; the `calibration_std` text states what σ_CAL
  covers (#9).

### Changed
- **[product]** The defaults are v0.5 (T49.1); the golden is v0.5 (T49.3,
  `docs/golden/MANIFEST_0.5.0.sha256`; 0.3 archived in `test_golden/archive/`).
- The image no longer installs the `download` extra (it brought dask and zarr,
  forbidden by `budget.yaml`); stage inputs from a separate install.
- `scripts/run_validation.sh` passes the golden's static ancillaries (frame
  table, area and event databases) when it has them.
- Venti pinned to the mgovorcin/venti-dev v0.5 stack (#3, #14, #15).
- Core install drops `opera-utils[disp]`; dependency and memory budget in
  `budget.yaml` (#5); peak memory v0.5 7.27 → 4.1 GB, gamma 6.0 → 3.8 GB (#8).
- UNR downloads go through geepers' verified downloader (#14).
- CI runs pre-commit on pull requests (#12).

### Fixed
- Frame-table overrides: Venti's zero-padded frame keys (`"08882"`) are found;
  before, no frame-table setting reached a run (#16).

### Golden pair: v0.5 vs 0.3 (the v0.5 golden, built and validated in `cal-disp:0.5.0-rc`)

v0.5 parameters + Venti frame table + area/event databases, as production will
run them:

| layer | median \|Δ\| | p99 \|Δ\| | max \|Δ\| | NaN share |
|---|---|---|---|---|
| `calibration` | 12.2 mm | 56.3 mm | 112.6 mm | 33.9% (water; unchanged) |
| `calibration_std` | 0.49 mm | 5.7 mm | 15.2 mm | 0% → 33.9%: NaN exactly where `calibration` is |
| `auxiliary/*` | 0 | 0 | 0 | identical |

The gamma file still reproduces the 0.3 golden in the image: every data layer
matches at 1e-6. Only the recorded parameter list differs, because the 0.3
record lacks the schema-v2 groups and is now read with v0.5 defaults.

### Open before the release
- ~~`calibration_std` NaN where the fit coverage is low~~ **decided
  2026-10-14:** σ_CAL follows the calibration (`uncertainty.extrapolate_sigma`
  in the v0.5 set; Venti #38): the nearest fitted value, inflated up to 3× with
  the distance, where the fit had too little weight. Golden pair: NaN masks of
  `calibration` and `calibration_std` equal, calibration and existing σ
  unchanged. To refine later (Venti todo).
- Per-frame `uncertainty.k_grid` (TS-G1) and tropo modes reach runs only
  through the frame table; the T46 runs did not pass it (k_grid 3.9 everywhere).
- Frame-table `frame.plate` is recorded but not used by any calibration step.

## [0.3.0] — gamma (golden 2026-10-06, not tagged)

Decision 0001 (`docs/decisions/0001-gamma-config.md`): 600 km window, tropo off,
unwrap correction off; golden rebuilt with `scripts/build_golden_output.sh`.
`configs/algorithm_parameters.yaml` and the model defaults still carry the
pre-decision values (30 km window, tropo on); runs that follow the README
therefore do not reproduce gamma. Use `configs/algorithm_parameters_gamma.yaml`.
