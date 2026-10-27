# 0001. Gamma 0.3 configuration

- **Status:** accepted 2026-10-06 (Venti plan T01.2)
- **Date:** 2026-10-06
- **Deciders:** OPERA ADT (Marin Govorcin), DISP-CAL lead (Talib Oliver-Cabrera)

## Why a decision is needed

The gamma release cannot be tagged until one configuration is agreed. Three
are in play, and none of the delivery documents states the algorithm
parameters (they live only in `configs/algorithm_parameters.yaml`):

| | window | tropo | unwrap correction | where it comes from |
|---|---|---|---|---|
| A. delivered gamma 0.3 | 600 km | on | **on** | delivery tar (cabrera) |
| B. current golden, upstream PR #21 | 600 km | on | off | `gamma-release` |
| C. Talib's validated stack | 30 km | off | off | F08882 stack; also GNSS-weighted, no global outlier mask |

## Evidence: the golden pair at real GNSS stations

F08882, 2022-01-11 → 2022-07-22. Calibrated DISP minus independent UNR
stations in LOS, mm (83 stations, actual motion between the two dates; MIDAS
rates give the same ranking). Same code (`gamma-release`) with one parameter
changed per row, except C. Sources: `trade studies/surface_fitting/
gnss_pair_validation/REPORT.md`, Venti `docs/benchmarks.md`.

| Configuration | bias | RMSE |
|---|---|---|
| A. delivered (600 km, tropo on, unwrap on) | −46.0 | 58.7 |
| B. current golden (600 km, tropo on, unwrap off) | +22.7 | 36.9 |
| **D. 600 km, tropo off, unwrap off** | **+1.1** | **21.4** |
| 30 km, tropo off, unwrap off | +0.4 | 19.6 |
| C. Talib's stack surface (30 km, tropo off, GNSS-weighted, no global mask) | +2.2 | 14.8 |
| *v0.5 candidate, for context (local-linear 50 km, tropo off)* | *+3.3* | *13.6* |

What drives it:

- **Unwrap correction on biases the product.** Venti's gamma corrector
  segments on the valid mask, shifts 5572 of 6733 regions and fits the surface
  to the shifted copy; calibrated DISP lands −9.3 cm from the GNSS grid
  (+0.5 cm with it off). It also never adds the shifts back into CAL.
- **The HRES troposphere correction over-corrects this frame.** Raw DISP −
  grid correlates +0.70 with the modelled differential delay; after the
  correction the residual correlates −0.37. Regressing the data on the model
  gives k ≈ 0.15, the product applies k = 1. Turning tropo off is worth
  15 mm RMSE here. Across four frames (trade study `tropo_check`) the full
  correction costs double-difference sill on every frame with < 0.4 km relief
  and helps only on Los Angeles (1.55 km relief).
- The window (600 vs 30 km) is worth under 2 mm once tropo is off.

## Recommendation

**D: 600 km, tropo off, unwrap off.** It is one flag away from the frozen
gamma algorithm (`apply_tropo_correction: false`), removes the two effects
that bias the product, halves the station RMSE of the current golden
(36.9 → 21.4 mm) and has the smallest bias of the gamma options. The 30 km
window and Talib's GNSS weighting are larger behaviour changes for under
2 mm more and belong to v0.5, where the local-linear surface already reaches
13.6 mm.

If tropo must stay on for gamma (for example because the delivery documents
promise it), B is the fallback; A should not ship.

## Consequences

- The golden changes (value change): rebuild it inside the Docker image
  (plan T03), re-run validation at 1e-6, update PR #21 and its description.
- The delivery documents describe the tropo layer as applied; the Product
  Spec §2 wording changes with the H items of `DOCUMENT_FIXES.md` (T04).
- v0.5 keeps tropo per frame through the frame-parameter table (off for
  F08882, F08886, F08622; stratified for F16940).

## Questions for Talib (with the proposed answer)

1. **Configuration:** A, B, C or D? *Proposed: D.*
2. **Unwrap correction off for gamma?** *Proposed: yes* (−9.3 cm bias with it
   on; the corrected estimator is in Venti, gated until TS-U1).
3. **Product version in the golden filename:** `v0.1` (delivered) or `v1.0`
   (ours)? *No proposal; needs the OPERA PGE convention.*
4. **Who updates the delivery documents** (memory 18 → 6.4 GB, Spec §2
   wording, URLs r2/r3)? *Proposed: ADT drafts the changes (plan T04) and
   Talib reviews.*

## Answers (2026-10-06)

The DISP-CAL lead (Talib Oliver-Cabrera) accepted all proposed changes:

1. **Configuration: D**, 600 km window, tropo off, unwrap off. Frozen in
   `scripts/build_golden_output.sh` (`apply_tropo_correction: false`); the
   golden is rebuilt inside the Docker image (Venti plan T03).
2. **Unwrap correction off** for gamma: yes.
3. **Product version in the golden filename:** no proposal was made; the
   golden keeps `v1.0`, what the code writes today. To confirm with the
   OPERA PGE team together with the product URLs (Venti plan T05).
4. **Delivery documents:** ADT drafts the changes (Venti plan T04), the
   DISP-CAL lead reviews.
