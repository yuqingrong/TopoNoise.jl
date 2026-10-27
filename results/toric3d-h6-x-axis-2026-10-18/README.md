# H+6 growth along the physical x axis

This study compares physical-x sequential preparation with the completed [diagonal H+6 study](../toric3d-h6-matched-2026-10-17/README.md). The physical direction is independent of the X/Z error channel: both error channels are included.

Status: complete. All 205 points reached 1,000 shots each (205,000 shared records and 615,000 protocol outcomes). The run finished on 2026-10-18 at 15:22:57 Asia/Shanghai, after 92.4 minutes. The final comparison is in [h6-direction-comparison.pdf](figures/h6-direction-comparison.pdf). This exploratory scan does not report a validated x-axis threshold.

All three feedback arms use the same immutable fault record at each point. The selected sizes are T=Ly=Lz=L=3,5,7,9,11, with 1,000 shots per point and the same physical-error-rate grids as the diagonal study. Timing and verification records are excluded from production. Directions use independent random streams.

| Property | Diagonal growth | Physical-x growth |
|---|---|---|
| Transfer per growth cell | H+6 CNOTs | H+6 CNOTs |
| Noisy transfer outputs | 13 L³ | 13 L³ |
| Initial toric slice | Honeycomb, 3 L² edges | Square yz plane, 2 L² edges |
| Logical qubits | 2 | 2 |
| Physical qubits | 3(L+1)L² | (3L+2)L² |
| Z logical distance | 2L | L |
| X logical distance | (L+1)L | (L+1)L |

The x-axis circuit first copies each old yz edge to its fresh counterpart, then applies H to each fresh x link and four CNOTs from that link to its incident fresh yz edges. All targets are fresh. Each completed prefix is the corresponding encoded cubic slab.

Feedback arms are final-only, one-step, and complete two-slice recovery. X one-step uses matching on square-slice plaquette defects and propagates only endpoint-edge corrections to the fresh plane. Z one-step pushes each retiring-star defect along its unique forward x edge; this differs from the diagonal matching-based endpoint rule. Two-slice recovery uses every physical check incident to the allowed two data slices and audits the full completed prefix. All arms share the sector's final decoder.

X starts with an ideal logical |00> seed. Z uses logical |++> seed-H noise followed by the same complete seed syndrome repair in all arms; logical seed errors and every trial are retained. The two-slice X observable is a preparation outcome fixed after the first corrected window, not a memory threshold in growth depth.

Equal L matches growth-cell count and transfer exposure, but does not match the finite code distance or seed geometry. This is an explicitly specified orientation comparison, not a direction-only control on identical finite codes. The initial scan is exploratory: diagonal fit windows are not assumed to bracket an x-axis threshold.

The bundle keeps consolidated data in `data`, a combined PDF in `figures`, and reproducible source, checkpoints, logs, and excluded pilots in `archive`. The existing diagonal study and manuscript are unchanged.

Verification passed before the source and sampling plan were frozen: 376 geometry assertions, 15,740 feedback assertions, 127 runner assertions, and 23 Python tests. All 18 timing-pilot points are excluded. Production contains 205 points, 205,000 shared records, and 615,000 protocol outcomes. Interrupted runs resume from checksummed batches.
