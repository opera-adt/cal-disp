# Golden dataset, v0.5.8

Frame F08882 (Houston), pair 2022-01-11 → 2022-07-22. The data live outside
git (`test_golden/`, 4.9 GB); `MANIFEST_0.5.8.sha256` pins every config,
input, static-ancillary and golden-output file
(`cd test_golden && sha256sum -c ../docs/golden/MANIFEST_0.5.8.sha256`).

| Item | Value |
|---|---|
| Configuration | v0.5 (`configs/algorithm_parameters_v05.yaml`, the defaults), with the Venti frame table (0.9-draft; for F08882: regions + stable-pixel fit mask) and the area (1.1-draft) and event (1.0-draft) databases in `input_data/static_ancillary/` |
| Built | 2026-11-01, inside `cal-disp:0.5.8-rc` (`sha256:19fa6488f3dc004d78cb2366f4820a99dbedabc29b62c65f3d4b78e1eaa69307`, version `0.5.8rc1`), with `scripts/build_golden_output.sh --skip-download`, in a directory named `golden_v058`, so the runconfig stored in the product names that path |
| Golden product | `OPERA_L4_DISP-CAL-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0_20261101T205953Z.nc` |
| Change vs 0.5.7 | `calibration` median +6.91 mm, p95 \|Δ\| 20.63, max 57.61; `calibration_std` median −0.19 mm, p95 \|Δ\| 1.58, max 10.71 (Venti `fc8c8f1`, venti-dev #55 and #57); other data layers identical |
| Validation in the image | passed at 1e-6, peak 5.77 GB (`test_golden/validation_0.5.8rc1.log`) |
| Change vs 0.3 | `calibration` median \|Δ\| 12.2 mm (max 113 mm); `calibration_std` NaN exactly where `calibration` is (CHANGELOG) |
| Previous golden (v0.5.7, 2026-10-28) | archived in `test_golden/archive/golden_v0.5.7_20261028/`; manifest `MANIFEST_0.5.7.sha256` |
| v0.5.5/0.5.6 golden (2026-10-22) | archived in `test_golden/archive/golden_v0.5.5_20261022/`; manifest `MANIFEST_0.5.5.sha256` |
| v0.5.4 golden (2026-10-22) | archived in `test_golden/archive/golden_v0.5.4_20261022/`; manifest `MANIFEST_0.5.4.sha256` |
| v0.5.2/0.5.3 golden (2026-10-20) | archived in `test_golden/archive/golden_v0.5.2_20261020/`; manifest `MANIFEST_0.5.2.sha256` |
| v0.5.1 golden (2026-10-16) | archived in `test_golden/archive/golden_v0.5.1_20261016/`; manifest `MANIFEST_0.5.1.sha256` |
| v0.5.0 golden (2026-10-15) | archived in `test_golden/archive/golden_v0.5.0_20261015/`; manifest `MANIFEST_0.5.0.sha256` |
| Gamma 0.3 golden (2026-10-06) | archived in `test_golden/archive/golden_v0.3_20261006/`; manifest `MANIFEST_0.3.0.sha256`; reproduced by `configs/algorithm_parameters_gamma.yaml` (data layers at 1e-6) |

Rebuilding: see the header of `scripts/build_golden_output.sh`. A change to
the configuration or to anything that moves the product is a new golden: rebuild
in the image, validate, regenerate this manifest.
