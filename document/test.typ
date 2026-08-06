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
    Complete verification report
  ]
]

#v(10pt)

= Overview

The complete suite checks the doubled-edge tensor and PEPS construction,
six-qubit local gates, graph and sequential-circuit representations, literal
Yao measurement/reset trajectories, an exact scalable sampler, raw
percolation observables, equal-weight marginalized component spins, Binder
cumulants and finite-size crossings, runnable command-line workflows, and all
generated figures. The verified run contains *1008 assertions* in *25 test
groups*.

#text(size: 8.3pt)[
  #table(
    columns: (2.55fr, 0.45fr, 2.55fr, 0.45fr),
    inset: 4pt,
    stroke: 0.5pt + rgb("b5b5b5"),
    fill: (_, row) => if row == 0 { rgb("e8eef7") },
    table.header(
      [*Test group*], [*Count*], [*Test group*], [*Count*],
    ),
    [Local toric-code tensor], [262],
    [Literal Yao measured trajectories], [8],
    [Finite PEPS structure], [25],
    [Scalable measured trajectories], [10],
    [Complete small-patch wavefunctions], [9],
    [Trajectory mismatch observables], [29],
    [$3 times 3$ patch], [6],
    [Marginalized component-spin observables], [25],
    [Unitary completion], [10],
    [Welford streaming accumulator], [3],
    [Toric-code local unitary gate], [19],
    [Streaming trajectory scans], [152],
    [Toric-code PEPS graph], [49],
    [Finite-size spanning crossings], [26],
    [Toric-code sequential circuit schedule], [122],
    [Marginalized Binder crossings], [24],
    [Executable Yao circuit], [26],
    [Runnable toric-code circuit example], [36],
    [Yao local-gate convention and small-patch states], [11],
    [Trajectory critical-scan CLI], [44],
    [Exact small-patch trajectory distributions], [25],
    [Trajectory scan figure], [21],
    [Virtual-bond trajectory model], [29],
    [Native Yao circuit figures], [13],
    [Toric-code vector figures], [24], [], [],
    [*Total*], [*1008*], [], [],
  )
]

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

= Graph and sequential-circuit tests

The `Toric-code PEPS graph` group verifies that the tensor network is keyed by
site coordinates, contains one graph edge for every shared internal virtual
index, and retains all four physical legs per site. It checks the expected
degrees of corners, edges, and bulk sites, confirms graph and PEPS dimensions
agree, and rejects invalid network input.

The sequential circuit uses one horizontal carrier per row and one vertical
carrier per column. Physical qubits are placed first in row-major site order
and `(E,N,W,S)` order within a site. The total qubit count is

$
n_("qubit") = 4 R C + R + C.
$

Sites are grouped into $R+C-1$ diagonals beginning at the lower-left corner.
The schedule tests check every site exactly once, the layer formula, gate and
wire metadata, deterministic ordering within a layer, and disjoint wires for
gates that may run in parallel. They also verify the analytic carrier
postselection probability

$
P_("post") = 2^(-(R+C)).
$

The executable Yao tests prepare carrier wires with Hadamards, apply the same
six-qubit local blocks as the sequential model, and compare small-patch states
against the tensor construction. Qubit counts, block structure, input
validation, scalar conventions, and the local gate's qubit ordering are all
checked without materializing a full circuit matrix.

= Virtual-bond error model

For an $R times C$ open patch, one independent Bernoulli-$p$ bit is sampled on
every internal and dangling virtual bond:

$
N_("bond") = R(C-1) + (R-1)C + 2R + 2C.
$

The `VirtualBondErrors` tests check all six storage groups independently:

#table(
  columns: (2.1fr, 1.35fr, 3fr),
  inset: 5pt,
  stroke: 0.5pt + rgb("b5b5b5"),
  fill: (_, row) => if row == 0 { rgb("e8eef7") },
  table.header([*Field*], [*Shape*], [*Location*]),
  [`horizontal_internal`], [$R times (C-1)$],
  [Between $(r,c)$ and $(r,c+1)$.],
  [`vertical_internal`], [$(R-1) times C$],
  [Between lower $(r+1,c)$ and upper $(r,c)$.],
  [`west_boundary`], [$R$], [Input end of each row carrier.],
  [`east_boundary`], [$R$], [Output end of each row carrier.],
  [`south_boundary`], [$C$], [Input end of each column carrier.],
  [`north_boundary`], [$C$], [Output end of each column carrier.],
)

The tests verify the exact total count, every expected shape, and each of the
six shape errors separately. At $p=0$ no bit is set, while at $p=1$ every bit
is set. A seeded statistical check selects one bit from each group over 5,000
draws: every marginal mean must lie within 0.025 of $p=0.31$, and every
pairwise joint mean must lie within 0.025 of $p^2$. Invalid, nonfinite error
rates are rejected.

= Measured trajectory tests

== Literal Yao measurement and reset

Horizontal carriers move west to east and vertical carriers move south to
north. The literal reference implementation uses four reusable physical
ancillas and $R+C$ carriers. The tested circuit timeline is:

+ Prepare every carrier in $|+ chevron.r$.
+ Apply west and south boundary $X$ events after preparation.
+ Reset the four physical ancillas to $|0000 chevron.r$.
+ Apply the six-qubit gate, measure physical `(E,N,W,S)` in the computational
  basis, record the result, and reset the same four ancillas.
+ Apply an eligible internal $X$ after the upstream measurement and before the
  downstream local gate.
+ Apply east or north boundary $X$ immediately after the final site on that
  carrier, before the analytic $chevron.l +|$ effect.

The literal trajectory tests verify the returned type, stored error object,
measurement-array shape, even local parity, and exact horizontal and vertical
mismatches. They also reject nonpositive `max_qubits` and a state-vector
allocation above the default 24-qubit guard.

Boundary events remain explicit even though

$
X |+ chevron.r = |+ chevron.r,
quad chevron.l +| X = chevron.l +|.
$

This identity explains why dangling errors are sampled metadata but do not
change the physical record or analytic postselection probability.

== Exact scalable sampler

The scalable sampler draws the same Born distribution without allocating a
state vector. At site $(r,c)$ it reads the incoming carrier bits
$(gamma,delta)$, samples $alpha$ uniformly, and computes

$
beta = alpha xor gamma xor delta,
quad (i,j,k,l) = (alpha,beta,gamma,delta).
$

The tests check local even parity and confirm that internal errors are applied
only after the upstream record. With fixed random seeds and fixed internal
errors, toggling all boundary metadata leaves the complete physical record
unchanged. Malformed error shapes are rejected.

== Exact small-patch distribution comparison

For $1 times 1$, $1 times 2$, and $2 times 1$, the test suite enumerates all
initial carrier bits and all sampled $alpha$ bits to construct the exact
classical record distribution. Independently, it applies the corresponding
Yao gates and marginalizes the carrier wires. The probability dictionaries
have identical support and agree to numerical tolerance, and both normalize
to one.

The comparison includes a fixed occupied internal bond for each nonsingleton
orientation. It then toggles every individual west, east, south, and north
dangling location on all three geometries. Every resulting Yao distribution
equals its boundary-free baseline.

= Mismatch and raw percolation observables

The doubled-edge record makes every internal virtual error directly visible.
With the array convention used by the implementation,

$
m^h_(r,c) = E_(r,c) xor W_(r,c+1),
quad
m^v_(r,c) = N_(r+1,c) xor S_(r,c).
$

The tests require these complete bit matrices to equal the stored horizontal
and vertical internal-error maps exactly. `trajectory_observables` refuses a
trajectory whose records and stored internal errors disagree.

The following raw quantities are checked:

- total sampled-error density over internal and boundary bonds;
- boundary-error density over the $2R+2C$ dangling legs;
- observable mismatch density over internal bonds;
- odd-plaquette frustration density;
- largest occupied-internal-bond cluster fraction from union-find;
- horizontal and vertical internal-error spanning flags.

A single occupied horizontal bond on a $2 times 2$ patch gives the expected
density, one frustrated plaquette, a two-site largest cluster, horizontal
spanning, and no vertical spanning. Boundary-only errors change the sampled
and boundary densities but leave all internal mismatch observables empty. At
$p=1$, every internal bond is occupied, the largest cluster fills the patch,
both directions span, and plaquette frustration is zero. The singleton case
has no internal mismatch or plaquette but one isolated-site cluster.

= Marginalized component-spin observables

For the spin analysis, the tests deliberately reinterpret every mismatched
internal doubled edge (`01` or `10`) as an erased relation. A matched edge
(`00` or `11`) requires its endpoint spins to agree. Union--find partitions
the $N=R C$ sites into matched components of sizes $n_a$. Open boundaries give
each component an independent sign $eta_a in {-1,+1}$, so

$
M = 1/N sum_a n_a eta_a.
$

Equal $1/2$ weights for both latent values of an erased edge make that edge
impose no spin constraint. The conditional even moments are therefore
calculated exactly without enumerating $2^(N_("component"))$ sign sectors:

$
chevron.l M^2 chevron.r_t = (sum_a n_a^2) / N^2,
quad
chevron.l M^4 chevron.r_t =
  (3 (sum_a n_a^2)^2 - 2 sum_a n_a^4) / N^4.
$

The tests exhaustively enumerate every component-sign sector for $1 times 1$,
$1 times 2$, $2 times 1$, and selected $2 times 2$ matched/erased patterns,
including a $3+1$ component decomposition. The enumerated moments are checked
against the analytic formulas, seeded component-sign samples converge to the expected
$chevron.l |M| chevron.r_t$, and repeated seeds reproduce the same estimate.
Changing stored latent-error metadata without changing the physical record
leaves all marginalized spin observables invariant. Nonpositive component-
sign sample counts are rejected.

For one component, both moments and $|M|$ equal one. For $N$ isolated sites,

$
chevron.l M^2 chevron.r = 1/N,
quad
chevron.l M^4 chevron.r = (3N^2-2N)/N^4,
quad
U_4 = 2/(3N).
$

= Streaming scans and finite-size crossings

== Online statistics

The Welford accumulator is tested on the sequence $(1,2,3,4)$, giving count
four, mean $2.5$, and the expected standard error. Finite-size scans retain
the seven raw online mean/standard-error pairs, marginalized $|M|$, $M^2$,
and $M^4$, horizontal-spanning batches, and paired spin-moment batches rather
than complete trajectories.

Seeded scans on sizes 2 and 3 check all output fields and reproducibility. At
$p=0$, injected, boundary, mismatch, frustration, and spanning observables
vanish; the largest cluster is one isolated site. At $p=1$, the three error
densities, largest-cluster fraction, and both spanning probabilities equal
one, while frustration remains zero. The spin sector gives $U_4=2/3$ at
$p=0$ and $U_4=2/(3L^2)$ at $p=1$. A derived RNG stream ensures that changing
`spin_samples` does not change any raw physical trajectory statistic. Binder
standard errors use a weighted delete-one-batch jackknife and are missing when
only one batch exists. A nonzero hand-calculated jackknife with unequal batch
counts verifies the weighting, and a mixed-component example distinguishes
the aggregate-moment Binder ratio from an average of trajectory ratios.
The same aggregate identity is asserted on a nontrivial generated scan point.
Invalid sizes, grids, shot counts, batch counts, and
component-sign sample counts are rejected.

== Spanning crossing classification and bootstrap

Adjacent-size horizontal-spanning curves are smoothed by pool-adjacent-
violators isotonic regression. An isolated sign-changing equality is a valid
crossing and is linearly interpolated; batch resampling supplies its confidence
interval. The crossing tests distinguish all three public statuses:

- `:ok`: a unique bracketed crossing with a sufficiently high valid-bootstrap
  fraction;
- `:unbracketed`: no interior sign change, including an endpoint-only contact;
- `:unstable`: identical curves, a sign-changing equal plateau, or fewer than
  80 percent valid bootstrap estimates.

At least two batches per scan point are required, preventing a meaningless
zero-width interval from one resampled batch. Invalid bootstrap counts and
confidence levels are also rejected.

== Binder crossing classification and bootstrap

The Binder ratio is formed after averaging conditional moments across shots,
never by averaging per-trajectory ratios:

$
U_4(p,L) = 1 -
  chevron.l M^4 chevron.r / (3 chevron.l M^2 chevron.r^2).
$

Binder curves are isotonic-smoothed in the non-increasing direction. Each
bootstrap replicate resamples paired $M^2$ and $M^4$ batches, recomputes the
ratio, and then locates one unique interior adjacent-size crossing. Exact
synthetic curves verify the estimate and confidence interval. No crossing is
`:unbracketed`; identical curves, sign-changing equal plateaus, multiple
crossings, and too few valid bootstrap estimates are `:unstable`. Seeded
nonconstant unequal-size batches verify paired resampling, while a noisy
synthetic scan deterministically exercises the below-80-percent status. Missing
spin points, one-batch input, invalid replicate counts, and invalid confidence
levels are rejected.

= Examples, CLI, and figure artifacts

The runnable circuit example is exercised with default and explicit lattice
sizes. It prints the lattice shape, qubit count, and diagonal layer schedule,
then writes nonempty SVG, PDF, and PNG circuit figures. Invalid dimensions and
argument counts return a usage error.

The trajectory scan CLI is run in a temporary directory with a seeded small
scan. Its required error-grid options are `--p-min`, `--p-max`, and
`--p-step`; the tests also exercise size, shot, batch, component-sign sample,
bootstrap, seed, and output-directory overrides. Missing options, a zero or
nondividing step, unsorted sizes, one batch, zero component-sign samples, and
unknown arguments are rejected.

The successful run checks all six output artifacts:

#table(
  columns: (2.4fr, 3.6fr),
  inset: 5pt,
  stroke: 0.5pt + rgb("b5b5b5"),
  fill: (_, row) => if row == 0 { rgb("e8eef7") },
  table.header([*Artifact*], [*Verification*]),
  [`trajectory_scan.csv`], [Exact 25-column raw and marginalized-spin header.],
  [`trajectory_crossings.csv`], [Exact 8-column header and preserved status.],
  [`trajectory_binder_crossings.csv`], [Separate 8-column Binder crossings.],
  [`trajectory_scan.svg`], [Nonempty XML vector figure.],
  [`trajectory_scan.pdf`], [Nonempty file with a PDF signature.],
  [`trajectory_scan.png`], [Nonempty raster figure.],
)

The six-panel scan figure retains all four raw panels and adds marginalized
$|M|$ and Binder panels with crossing markers. A legacy scan without spin
points still renders its original four-panel layout. A negative-Binder
regression checks that the vertical range expands instead of clipping valid
data, and unavailable one-batch uncertainties do not render as zero-width
error bars. All three save formats
are checked. Native Yao figures are validated separately for $1 times 1$ and
$2 times 2$ circuits. CairoMakie PEPS and sequential-circuit figures are
checked for site, direction, carrier, gate, state, and layer labels before SVG
and PDF export.

#block(fill: rgb("f5f1e8"), inset: 8pt, radius: 2pt)[
  *Scope.* The spin analysis marginalizes `01/10` edges with equal weights,
  making them deleted constraints between independently oriented matched
  components. Its Binder cumulant describes this bond-diluted component-spin
  ensemble, not a decoder or finite-temperature random-bond Ising model.
]

= Core PEPS test helpers

The original tensor-network tests use five focused helper routines:

#text(size: 8.3pt)[
  #table(
    columns: (1.7fr, 4fr),
    inset: 2.5pt,
    stroke: 0.5pt + rgb("b5b5b5"),
    fill: (_, row) => if row == 0 { rgb("e8eef7") },
    table.header([*Helper*], [*Role*]),
    [`ordered_physical_indices`], [Returns physical indices in site and east/north/west/south order.],
    [`contract_state`], [Fully contracts a small PEPS while retaining all physical indices.],
    [`is_allowed_configuration`], [Checks local parity and duplicated-bond agreement.],
    [`projected_amplitude`], [Projects a basis configuration and contracts only the virtual network.],
    [`double_layer_norm2`], [Contracts the exact bra-ket norm without constructing the dense wavefunction.],
  )
]

#text(size: 8.3pt)[*Running the tests.* Run
`julia --project=. -e 'using Pkg; Pkg.test()'` with Julia
1.12.6, ITensors 0.9.30, ITensorNetworks 0.21.5, Yao 0.9.3, and CairoMakie
0.15.13: all *1008 assertions* passed. The Typst build and `git diff --check`
passed; independent review found no critical or important issues; and
`document/note.tex` was unchanged.]
