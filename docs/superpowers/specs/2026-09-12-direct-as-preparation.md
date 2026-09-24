# Direct As logical-state preparation

The CSS convention remains `A_s = product(Z)` and `B_p = product(X)`.
The goal is to prepare the selected logical state using the existing As
plaquette CNOTs once, removing the As inverse followed by a Bp reconstruction.

## Circuit and derivation

Let `R` be the fresh representatives of the existing As blocks and `F` the
remaining input qubits. Every representative is untouched by preceding core
blocks. Its block H can therefore commute backwards and cancel the H from
the initial all-plus preparation. The natural As preparation becomes

```math
|+_L\rangle = C\bigl(|0\rangle_R\otimes|+\rangle_F\bigr),
```

where `C` is the original ordered sequence of As CNOTs, with no block H gates.
The constructor checks the fresh-representative invariant before cancelling H.

Conjugate the canonical logical Z through this CNOT core:

```math
C^\dagger\bar Z C=\prod_{q\in Q}Z_q,\qquad T=Q\cap F.
```

The factors on representatives already have eigenvalue +1. Preparing even
Z parity on `T` therefore fixes the output logical Z to +1. Choose a pivot
`t` in `T`, initialize `t` and all representatives in `|0>`, and apply H to
every qubit in `F` except `t`. Apply `CNOT(q,t)` for each other `q` in `T`,
in increasing qubit order, then run `C` once. The input is the normalized
even-parity projection of the natural input. Thus the output is

```math
|0_L\rangle\propto(I+\bar Z)\prod_s(I+A_s)|+\rangle^{\otimes n}.
```

This is a deterministic preparation from all-zero input; its unitary circuit
does not implement a nonunitary projector on arbitrary input states.

The pivot minimizes `(maximum Manhattan distance to T, total Manhattan
distance to T, qubit index)` lexicographically. This is deterministic and
keeps the distance-3 parity CNOTs local. At larger distances the parity
CNOTs can be nonlocal. No SWAP routing or ancillary qubits are inserted.

For `:one`, append the canonical logical X string to the zero preparation.
For `:plus`, apply H to every free input and omit the parity CNOTs. For
`:minus`, append the canonical logical Z string to the plus preparation.
Both boundary orientations use the same rule on their respective core and
canonical logical support. The Bp schedules are unchanged for all four states.

## Distance-3 example

For the default `:x_ns` orientation, qubits are numbered by rows:

```text
1 2 3
4 5 6
7 8 9
```

The representatives are `[4,5,7,9]`. The free-input parity support is `[1,2,3]`,
and the pivot is 2. Starting from nine zeros, execute:

```text
H(1), H(3), H(6), H(8)
CNOT(1,2), CNOT(3,2)
CNOT(1,4)
CNOT(2,5), CNOT(3,5), CNOT(6,5)
CNOT(4,7), CNOT(5,7), CNOT(8,7)
CNOT(6,9)
```

This is 14 gates, compared with 54 in the previous As zero preparation.

## Inspection and noise

The initial H gates remain in the `:logical` block. The parity CNOTs occupy
a separate `kind=:logical_sector`, `source_check=:logical_parity` block.
Its `source_support` is `T`, and its `representative` is the parity target.
As plaquette blocks preserve their check indices, supports, representatives,
execution order, and incoming CNOT directions.

Gate-layer noise includes both initialization blocks on their literal active
qubits. Plaquette noise includes only completed source-check blocks. Pauli
frame replay, Yao replay, and the decoder consume the stored circuit without
construction-specific changes. Earlier As noise results need a new run.

This preparation is not fault tolerant: for example, an X on qubit 2 just
after `CNOT(1,2)` propagates to `X2 X5 X7`, which has zero stabilizer syndrome
but anticommutes with logical Z. A Z at the same location propagates to
`Z2 Z3`, with a nonzero Bp syndrome. Cancelling gates reduces fault locations
but does not guarantee that an individual preparation fault is correctable.

## Verification requirements

- Preserve the distance-3 CNOT sequence and the direct 14-gate zero circuit.
- Verify all stabilizers and the signed logical observable for all four states,
  both orientations, and odd distances from 3 through 15.
- Independently compare the distance-3 statevector with the projector formula
  and the Bp preparation.
- Check the new preparation fault locations and hand-derived propagated
  Pauli supports against the Yao syndrome oracle.
- Run the complete package tests, including decoder, scans, and plotting.
