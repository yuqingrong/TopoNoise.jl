#set page(
  paper: "a4",
  margin: (x: 24mm, y: 22mm),
  numbering: "1",
)
#set text(font: "New Computer Modern", size: 10.5pt)
#set par(justify: true, leading: 0.68em)
#set heading(numbering: "1.1")
#show raw: set text(font: "DejaVu Sans Mono", size: 9pt)

#align(center)[
  #text(size: 18pt, weight: "bold")[TopoNoise.jl Tests]
  #v(4pt)
  #text(size: 10pt, fill: rgb("555555"))[
    Summary of `tests`
  ]
]

#v(10pt)

= Overview

The test suite validates the local doubled-edge toric-code tensor, the
connectivity and normalization of finite PEPS patches, complete wavefunctions
for small systems, a memory-safe exact check of a $3 times 3$ patch, and the
completion of the local isometry into a six-qubit unitary gate. The suite
contains *331 assertions* in six test groups.

#table(
  columns: (2.2fr, 0.7fr, 3.6fr),
  inset: 6pt,
  stroke: 0.5pt + rgb("b5b5b5"),
  fill: (_, row) => if row == 0 { rgb("e8eef7") },
  table.header(
    [*Test group*], [*Assertions*], [*Purpose*],
  ),
  [Local toric-code tensor], [262],
  [Check every local tensor entry, the virtual matrix $W$, and the local isometry.],
  [Finite PEPS structure], [25],
  [Check dimensions, indices, bonds, input validation, and scalar types.],
  [Complete small-patch wavefunctions], [9],
  [Enumerate all amplitudes for $1 times 1$, $1 times 2$, and $2 times 2$.],
  [$3 times 3$ patch], [6],
  [Check selected amplitudes and the exact double-layer norm without forming the dense state.],
  [Unitary completion], [10],
  [Check generic nullspace completion, two-sided unitarity, and invalid inputs.],
  [Toric-code local unitary gate], [19],
  [Check the rank-12 gate, scalar types, unitarity, and the exact $|0000 chevron.r$ input slice.],
)

= Local tensor tests

The local tensor uses physical indices $(i,j,k,l) = (E,N,W,S)$ followed by
virtual indices $(alpha,beta,gamma,delta) = (E,N,W,S)$. Binary values 0 and 1
occupy Julia positions 1 and 2. Its normalized definition is

$
T^(i j k l)_(alpha beta gamma delta)
= cases(
  1 / sqrt(2) delta_(i alpha) delta_(j beta) delta_(k gamma) delta_(l delta)
    & alpha + beta + gamma + delta = 0 quad "mod" quad 2,
  0 & "otherwise",
)
$

The tests verify:

- the array has eight indices, each of dimension 2;
- exactly eight entries are nonzero;
- all $2^8 = 256$ entries equal $1 / sqrt(2)$ exactly when the copy constraints
  and even-parity rule hold, and equal zero otherwise;
- integer storage is rejected because it cannot represent the normalization;
- requesting `ComplexF64` produces a complex-valued tensor;
- contracting all physical indices produces the expected virtual parity matrix;
- the tensor is an isometry from $(alpha,beta)$ to
  $(gamma,delta,i,j,k,l)$.

For the last check, each physical leg is contracted with the uniform covector
$(1,1)$. The copy constraints transfer the physical values to their
corresponding virtual values:

$
tilde(W)_(alpha beta gamma delta)
= sum_(i,j,k,l) T^(i j k l)_(alpha beta gamma delta).
$

The first two and last two virtual legs are fused into row and column indices,
respectively. Because $T$ is already normalized, the contraction directly gives

$
W = 1 / sqrt(2) mat(
  1, 0, 0, 1;
  0, 1, 1, 0;
  0, 1, 1, 0;
  1, 0, 0, 1;
).
$

For the isometry test, the stored order
$(i,j,k,l,alpha,beta,gamma,delta)$ is permuted to
$(gamma,delta,i,j,k,l,alpha,beta)$ and reshaped into a $64 times 4$ matrix

$
T_((gamma delta i j k l),(alpha beta)).
$

Each $(alpha,beta)$ input column contains two nonzero entries of magnitude
$1 / sqrt(2)$, while distinct columns have disjoint support. Therefore

$
T^dagger T = I_4.
$

= Unitary completion and local gate tests

== Generic unitary completion

The generic routine `unitary_completion` is tested with the fixed complex
$4 times 2$ isometry

$
A = mat(
  1 / sqrt(2), 0;
  i / sqrt(2), 0;
  0, 1;
  0, 0;
).
$

The implementation computes an orthonormal basis $N$ for
$ker(A^dagger)$ and returns the square matrix $U = (A quad N)$. The tests
check that `size(U) == (4, 4)`, its first two columns equal $A$ exactly, both
$U^dagger U = I_4$ and $U U^dagger = I_4$ hold numerically, and a square
unitary input is returned unchanged. Separate error assertions reject a wide
matrix, integer storage, non-orthonormal columns, and negative absolute or
relative tolerances.

Only the leading isometry columns are fixed by the construction. The remaining
orthogonal basis is non-unique, so the tests check its mathematical properties
rather than pinning particular nullspace entries.

== Toric-code local unitary gate

For the gate construction, the stored local tensor is directly reshaped using
the diagram orientation

$
T_((i j k l alpha beta),(gamma delta)) in CC^(64 times 4),
quad T^dagger T = I_4.
$

This differs from the $(alpha,beta)$ input cut checked in the preceding local
tensor section. The copy and parity symmetry of this tensor makes both cuts
isometric. Completing the diagram-oriented matrix produces a $64 times 64$
unitary, returned as a rank-12 array with output indices followed by input
indices:

$
(i,j,k,l,alpha,beta;
  "phy"_i,"phy"_j,"phy"_k,"phy"_l,gamma,delta).
$

The gate tests run for both `Float64` and `ComplexF64`. For each scalar type,
they verify all twelve dimensions are two, the requested element type is
preserved, and the reshaped gate satisfies both

$
U^dagger U = I_64, quad U U^dagger = I_64.
$

#block(breakable: false)[
  Julia position 1 represents qubit value zero. Fixing the four physical input
  indices to position 1 is therefore tested with exact array equality:

  $
  U_((i j k l alpha beta),(0 0 0 0 gamma delta))
  = T^(i j k l)_(alpha beta gamma delta).
  $
]

The same equality is checked separately for all four computational-basis
choices of $(gamma,delta)$. Finally, requesting an integer-valued gate is
verified to raise an error because the normalized tensor cannot be stored in
an integer array.

= Finite PEPS structure tests

The constructor `toric_code_peps(rows, cols)` is checked on a $3 times 3$
patch. The structural tests verify:

- zero rows or columns are rejected;
- integer scalar storage is rejected because the normalized state requires
  floating-point square roots;
- `size(peps)` and dimension-specific size queries return the lattice shape;
- a corner, edge, and bulk tensor have 6, 7, and 8 remaining indices,
  respectively: four physical indices plus the internal virtual bonds;
- all 36 physical qubit indices are unique;
- each horizontal or vertical neighboring pair shares exactly one virtual
  index;
- the $3 times 3$ patch has 12 unique internal bonds;
- non-neighboring corner tensors share no indices;
- `ComplexF64` construction preserves the requested tensor element type.

= Complete small-patch wavefunctions

The tensors are fully contracted for patches of size $1 times 1$, $1 times 2$,
and $2 times 2$. Let

$
V = R C, quad
E = R(C - 1) + (R - 1)C, quad
B = 2R + 2C.
$

Here $V$ is the number of vertices, $E$ is the number of internal bonds, and
$B$ is the number of virtual boundary legs terminated by $|+ chevron.r$ vectors.
There are

$
N_("allowed") = 2^(E + B - V)
$

allowed configurations, each with amplitude

$
psi_("allowed") = 2^(-(E + B - V)/2).
$

To see the count, assign one binary variable to every internal bond and every
boundary leg. This gives $E+B$ freely chosen bits before imposing the local
tensor constraints. Each of the $V$ vertices contributes one independent
even-parity equation over $bb(F)_2$. The equations are independent on an open
patch because every connected component reaches a boundary leg: in any
nonempty sum of vertex equations, a boundary leg incident on a selected
boundary vertex appears exactly once, so the sum cannot vanish. The allowed
bit space therefore has dimension $(E+B)-V$, and hence contains
$2^(E+B-V)$ configurations. Once a bond-bit assignment is fixed, the copy
constraints uniquely determine all four physical bits at every site, so no
additional multiplicity appears.

For each small patch, the complete amplitude table is checked to ensure that:

- the number of nonzero amplitudes equals $N_("allowed")$;
- an amplitude is nonzero exactly when every site has even parity and the two
  physical copies of every internal bond agree;
- every allowed configuration has amplitude $psi_("allowed")$;
- the sum of squared amplitudes is one.

= Exact $3 times 3$ patch tests

A dense $3 times 3$ wavefunction would contain $2^36$ amplitudes, so the test
does not construct it. Instead, individual physical configurations are
projected locally before the remaining virtual network is contracted.

For this patch,

$
V = 9, quad E = 12, quad B = 12,
quad E + B - V = 15,
$

so every allowed configuration must have amplitude

$
2^(-15/2) = 1 / sqrt(32768).
$

The following configurations are tested:

#table(
  columns: (1.7fr, 3.8fr, 1.3fr),
  inset: 6pt,
  stroke: 0.5pt + rgb("b5b5b5"),
  fill: (_, row) => if row == 0 { rgb("e8eef7") },
  table.header([*Configuration*], [*Reason*], [*Expected*]),
  [All-zero state], [Every vertex has even parity and all bond copies agree.],
  [$2^(-15/2)$],
  [Single plaquette loop], [Four internal bonds form a closed loop.],
  [$2^(-15/2)$],
  [Boundary-to-boundary string], [The open string terminates on the x-polarized virtual boundary.],
  [$2^(-15/2)$],
  [Odd vertex], [One vertex violates the local even-parity constraint.], [0],
  [Mismatched bond copies], [The two physical copies of one internal bond disagree.], [0],
)

Finally, the norm is evaluated exactly as a double-layer tensor network. Each
local bra tensor has primed virtual indices and unprimed physical indices, so
the physical legs contract locally while bra and ket virtual networks remain
distinct. Contracting the resulting $3 times 3$ double layer gives

$
chevron.l psi | psi chevron.r = 1.
$

= Test helper operations

The test file uses five small helper routines:

#table(
  columns: (1.7fr, 4fr),
  inset: 6pt,
  stroke: 0.5pt + rgb("b5b5b5"),
  fill: (_, row) => if row == 0 { rgb("e8eef7") },
  table.header([*Helper*], [*Role*]),
  [`ordered_physical_indices`], [Returns physical indices in site and east/north/west/south order.],
  [`contract_state`], [Fully contracts a small PEPS while retaining all physical indices.],
  [`is_allowed_configuration`], [Checks local parity and duplicated-bond agreement.],
  [`projected_amplitude`], [Projects a basis configuration and contracts only the virtual network.],
  [`double_layer_norm2`], [Contracts the exact bra-ket norm without constructing the dense wavefunction.],
)

= Running the tests

From the repository root, run
`julia --project=. -e 'using Pkg; Pkg.test()'`. The tested configuration is
Julia 1.12.6 with ITensors 0.9.30.
