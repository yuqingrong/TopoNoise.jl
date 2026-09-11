# Local As/Bp plaquette-growth encoder

## Goal

Replace the current RREF-derived rotated-planar encoder schedule with the
local sequential constructions shown in `document/image/2DcircuitBp.typ` and
`document/image/2DcircuitAs.typ`. The replacement must preserve the CSS
convention

\[
A_s=\prod Z,\qquad B_p=\prod X,
\]

and prepare the same selectable single-logical-qubit states. It must model
the intended circuit-level propagation: a newly introduced representative is
fresh at its block, while a diagonal same-family plaquette may reuse one
older corner qubit in its prescribed non-representative role.

The existing high-shot plots describe the RREF schedule only. They are not
results for this local growth construction and will not be reused after the
replacement.

## Authoritative d=3 circuits

The literal gate lists in the two Typ files are golden schedules. Site names
use the diagram's `siteXY` convention, which maps to the row-major code
coordinate `(row=d-Y+1, column=X)`. Under that map the diagrams are the
direct template for `boundary_orientation=:x_ew`; `:x_ns` applies the same
template after a 90-degree clockwise rotation. The Bp construction is:

| Representative | Gates after its H | Block kind |
|---|---|---|
| `site21` | `21 -> 11`, `21 -> 22`, `22 -> 12` | bulk |
| `site13` | `13 -> 12` | boundary, weight 2 |
| `site31` | `31 -> 32` | boundary, weight 2 |
| `site33` | `33 -> 32`, `33 -> 23`, `23 -> 22` | bulk |

The As construction is:

| Representative | Gates after its H | Block kind |
|---|---|---|
| `site21` | `11 -> 21` | boundary, weight 2 |
| `site22` | `12 -> 22`, `13 -> 22`, `23 -> 22` | bulk |
| `site33` | `23 -> 33` | boundary, weight 2 |
| `site31` | `21 -> 31`, `22 -> 31`, `32 -> 31` | bulk |

For its sequential d=3 execution, the As blocks are ordered `site21`,
`site22`, `site31`, `site33`: this is the southwest-to-northeast geometric
order and agrees with the time offsets in the Typ drawing. The table retains
the source-file grouping so every literal gate still has a visible home.

These schedules, including the order *within a plaquette block*, CNOT
direction, and fresh red representatives, supersede RREF pivots and reduced
supports. Plaquette blocks execute to completion in the stated deterministic
southwest-to-northeast order; gate heights in the Typ drawing are a visual
layout, not a second overlapping execution clock. The d=3 unit tests will
compare the stored elementary operations exactly with these lists.

## Geometry and growth rule

The patch remains an odd square rotated planar code with one logical qubit.
Same-family plaquettes are checkerboard-separated: no two share an edge, but
diagonal plaquettes may share one corner data qubit. Bulk checks have weight
four and exposed boundary checks have weight two.

For d=3, each construction follows the literal southwest-to-northeast Typ
order exactly, including after the orientation rotation. For d>=5, the
generator uses a deterministic reverse-shell order: it peels removable
checks by their still-unique data-qubit count and check index, then reverses
that list for execution. This preserves a fresh local representative and,
for a bulk check, a fresh adjacent continuation. `geometric_order` records
execution order; it is not required to sort geometric keys.

A path contains one representative and either one boundary edge or three
local CNOT edges forming a bulk tree. The first operation is H on the
representative.

- Bp paths use outgoing CNOT flow from the fresh representative, with the
  diagram's prescribed directed continuation through an older shared corner.
- As paths use incoming CNOT flow into the fresh representative.
- Every representative must be absent from all preceding plaquette-block
  operations of that construction. A shared diagonal corner may occur in a
  later block only as a non-representative.
- No algebraic row combination, pivot selection, nonlocal routing SWAP, or
  `reduced_support` is permitted. The block support is its actual local
  check support and its provenance is exactly that one source check.

The d=3 templates establish the parity/orientation choice for the general
enumerator. It must produce the same supports, representatives, CNOTs, and
operations as the tables above. The general rule is accepted only after
tableau verification at d=3, d=5, and d=7 for both constructions and both
boundary orientations.

## State and noise semantics

The Bp core maps the all-zero product state to \(|0_L\rangle\). The As
core, preceded by H on every data qubit, maps it to \(|+_L\rangle\). The
encoder retains both literal cores but applies a local state-specific basis
conversion when necessary: `As†; H_all; Bp` after an As core produces
\(|0_L\rangle\), and `Bp†; H_all; As` after a Bp core produces
\(|+_L\rangle\). Logical X or Z strings then produce \(|1_L\rangle\) or
\(|-_L\rangle\), respectively. Thus `construction` changes the circuit,
not the requested final logical state.

Circuit Pauli noise remains unchanged operationally: a `:plaquette` fault is
sampled after each completed local block on exactly that block's physical
support. Consequently, for Bp an older X error sitting on a later target is
not propagated by that later CNOT, while an allowed Z-on-target component
can propagate back toward its new control. The dual directional statement
holds for As. `:gate_layer` continues to use the exact same stored elementary
layers.

Final syndrome and decoding remain CSS algebra:

\[
s_{A_s}=x\text{-frame parity},\qquad
s_{B_p}=z\text{-frame parity},\qquad
\ell_X=x\cdot\bar Z,\qquad
\ell_Z=z\cdot\bar X.
\]

## Implementation boundaries

Modify the encoder subsystem and its tests. Geometry may be adjusted only
where required for the authoritative d=3 support tables and their rotated
orientation; the decoder, Pauli-frame propagation, scan API, comparison
workflow, and plotting API must consume the new canonical schedule without
special cases.

The current RREF helper functions and all metadata whose meaning relies on
row reduction (`source_rows` with multiple elements, `reduced_support`, and
algebraic pivots) will be removed or replaced with local equivalents. Public
block inspection remains available, now exposing the literal source support,
representative, local direction, and elementary layers.

## Verification

Develop test-first.

1. Add exact d=3 operation-sequence tests for both Typ schedules.
2. Check bulk/boundary support weights, same-family diagonal-only sharing,
   source support equality, local nearest-neighbour CNOTs, deterministic
   sweep order, and fresh representatives at d=3, 5, and 7.
3. Retain and rerun stabilizer commutation/rank/logical-string tests.
4. Verify each construction's tableau and exact d=3 Yao statevector against
   all A_s, B_p, and selected logical observables.
5. Retain frame-versus-Yao fault replay and final-syndrome tests, then add
   directed single-fault propagation cases demonstrating Bp-X and As-Z
   target suppression at the local-block boundary.
6. Run the full package test suite and `git diff --check`.
7. Only after the circuit tests pass, regenerate channel-comparison data;
   label old RREF-derived result directories as superseded rather than mixing
   them with new-circuit results.

## Non-goals

This change does not add repeated syndrome rounds, measurement noise,
circuit-aware matching, periodic patches, holes, more logical qubits, or a
threshold claim. It does not modify the definitions of A_s/Z, B_p/X,
logical-X, logical-Z, or the Pauli-frame decoder sectors.
