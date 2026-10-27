# Golden dataset, gamma 0.3

Frame F08882 (Houston), pair 2022-01-11 → 2022-07-22. The data live outside
git (`test_golden/`, 4.9 GB); `MANIFEST_0.3.0.sha256` pins every config,
input and golden-output file (`cd test_golden && sha256sum -c ../docs/golden/MANIFEST_0.3.0.sha256`).

| Item | Value |
|---|---|
| Configuration | 600 km window, tropo **off**, unwrap off, SET on ([decision 0001](../decisions/0001-gamma-config.md)) |
| Built | 2026-10-06, inside `cal-disp:0.3.0-rc2` (`sha256:95d5fffde592c24730ad5e717c4c39b71d99d3b91953025e8213e00536f3b3c0`), from `gamma-release` `ac9616f`, with `scripts/build_golden_output.sh --skip-download` |
| Golden product | `OPERA_L4_DISP-CAL-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0_20261006T224844Z.nc` |
| Validation in the image | passed at 1e-6, peak 5.65 GB (`test_golden/validation_0.3.0.log`) |
| Validation on aurora, pixi `dev` (Venti `2a7e61f`, GDAL 3.13.3) | passed at 1e-6 (`test_golden/validation_0.3.0_aurora_pixi.log`) |
| Previous golden (tropo on, 2026-09-30) | archived in `test_golden/archive/` |

Rebuilding: see the header of `scripts/build_golden_output.sh`. A change to
the configuration or to anything that moves the product is a new golden: rebuild
in the image, validate, regenerate this manifest.
