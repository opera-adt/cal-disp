# cal-disp

`cal-disp` calibrates OPERA DISP-S1 surface-displacement products to an
absolute reference. A DISP-S1 product gives the line-of-sight (LOS)
displacement between a reference and a secondary date, relative to an
arbitrary reference pixel. `cal-disp` fits a smooth *calibration surface*
that ties that relative field to the UNR gridded GNSS solution, and writes
it as the OPERA L4 **DISP-CAL-S1** product. Subtracting the calibration layer
from the DISP displacement gives absolute LOS displacement in the GNSS
reference frame (IGS20).

## What the workflow does

1. **Inputs** — one DISP-S1 product, the DISP-S1-STATIC line-of-sight and
   DEM layers for the frame, the UNR gridded GNSS time series (constant
   velocity or time-variable grid), and optionally OPERA TROPO zenith-delay
   products for both dates.
2. **GNSS reference** — the GNSS displacement between the two dates is
   projected onto the radar LOS and interpolated over the frame.
3. **Corrections** — tropospheric delay (from TROPO) and the solid-earth
   tide (from the DISP product) are removed before the fit and added back to
   the calibration surface, so the surface cannot absorb them.
4. **Surface fit** — [Venti](https://github.com/opera-adt/Venti) fits the
   calibration surface to `displacement − GNSS` on a downsampled grid and
   upsamples it back to the 30 m DISP posting.
5. **Product** — the surface, its uncertainty, identification and metadata
   groups, and a browse image are written as a compressed NetCDF product
   named like the input, e.g.
   `OPERA_L4_DISP-CAL-S1_IW_F08882_VV_20220111T002651Z_20220722T002657Z_v1.0_<production time>.nc`.

## Where to go next

- [Quickstart](quickstart.md) — install, stage inputs, configure, run,
  validate.
- [Workflow walkthrough](notebooks/00_calibration_walkthrough.ipynb) — the
  workflow executed step by step in a notebook, ending with a comparison
  against the golden dataset.
- [Fix inspection notebooks](notebooks/01_unwrap_cycle_length.ipynb) — one
  runnable notebook per defect found in the gamma-release review, showing
  the behaviour before and after the fix.
- [Development](development.md) — tests, linting, lock files, the golden
  dataset, and the release checklist.
- [API reference](api.md).

## Command line

```text
cal-disp config           Create a calibration workflow configuration file.
cal-disp download         Stage inputs: disp-s1, unr, tropo, burst-bounds.
cal-disp run              Run the calibration workflow for a runconfig.
cal-disp validate         Validate a DISP-CAL product against a reference.
cal-disp validate-golden  Re-run the calibration on golden inputs and compare.
```

`opera_cal-disp` is the same CLI under the name used in the delivery
instructions.
