# H+5 Z-boundary ambiguity diagnostics

This separate result set audits the current three-dimensional direct-projector H+5 upper-boundary instrument on cubic `L=3,5,7`. It does not alter or reinterpret any earlier production checkpoint, selection, figure, report, or LaTeX note.

The exhaustive final-layer audit compares the existing unit-length MWPM recovery with a circuit-aware finite-size control. The current MWPM has 34, 131, and 276 single-fault logical-Z failures at `L=3,5,7`. The best decoder using only the final star syndrome and observed upper parity still has 8, 14, and 20 failures. Those nonzero counts are an information-theoretic floor for the audited observable set, so the results are `diagnostic_only` and no threshold fit is attempted.

The exact upper outcome string does not reveal a fault location. The clean upper branch is uniform over even bit strings and the observed branch is its XOR with the fault-induced upper syndrome. Conditioning on the latter therefore leaves the observed string uniform within one parity class; only its parity carries fault information.

`parity-erasure-control.csv` reports a 1000-record paired diagnostic at `p=10^-4,10^-3` for all four protocols. Odd parity is treated as a physically available erasure flag, with rejected shots and acceptance intervals reported explicitly. This is a finite-size postselected control and cannot support a deterministic threshold claim.

The measured counts are:

| L | p | protocol | accepted even | rejected odd | baseline → control failures | accepted baseline/control failures |
|---:|---:|:---|---:|---:|:---|---:|
| 3 | 0.0001 | none / one_step | 980 | 20 | 14 → 10 | 0 |
| 3 | 0.0001 | window_complete / full_prefix | 992 | 8 | 6 → 2 | 0 |
| 3 | 0.001 | none | 817 | 183 | 109 → 87 | 7 |
| 3 | 0.001 | one_step | 817 | 183 | 101 → 79 | 7 |
| 3 | 0.001 | window_complete / full_prefix | 938 | 62 | 55 → 32 | 22 |
| 5 | 0.0001 | none | 921 | 79 | 55 → 50 | 2 |
| 5 | 0.0001 | one_step | 921 | 79 | 51 → 44 | 2 |
| 5 | 0.0001 | window_complete / full_prefix | 985 | 15 | 14 → 5 | 2 |
| 5 | 0.001 | none | 531 | 469 | 387 → 352 | 99 |
| 5 | 0.001 | one_step | 531 | 469 | 351 → 304 | 83 |
| 5 | 0.001 | window_complete / full_prefix | 830 | 170 | 160 → 85 | 52 |
| 7 | 0.0001 | none | 777 | 223 | 160 → 138 | 15 |
| 7 | 0.0001 | one_step | 777 | 223 | 155 → 127 | 14 |
| 7 | 0.0001 | window_complete / full_prefix | 947 | 53 | 49 → 13 | 11 |
| 7 | 0.001 | none | 486 | 514 | 593 → 592 | 246 |
| 7 | 0.001 | one_step | 486 | 514 | 569 → 554 | 240 |
| 7 | 0.001 | window_complete / full_prefix | 723 | 277 | 303 → 184 | 102 |

The accepted-failure count is identical for the baseline and control in every row: the circuit-aware correction acts on the odd-parity ambiguity classes, while even-parity postselection removes that first-order sector. The remaining accepted failures include higher-order events. Unknown multi-fault signatures use the unchanged MWPM recovery and are counted in the CSV's fallback column.

Logical Z here is a state-independent residual Pauli diagnostic. It is not `|00_L>` infidelity; that prepared-state failure is associated with residual logical X.

A deterministic positive threshold requires a later physical-boundary or dimensional redesign that supplies additional local information. Decoder tuning cannot reconstruct observables that the present instrument does not measure.

## Reproduction

Use the single-thread environment recorded in `runtime-environment.json`, then run `run_audit.jl`, `run_control.jl`, and `apply_guard.py` in that order. The Julia drivers write atomically after all assertions pass. The Python driver independently cross-checks JSON, group CSV, control CSV, and Markdown counts before writing `guard.json`.

The completed small-lattice scan remains at `../toric3d-h5-diagonal-feedback-2026-10-13/preliminary-l3-7`. The deliberately stopped full production checkpoint remains preserved and was not resumed.
