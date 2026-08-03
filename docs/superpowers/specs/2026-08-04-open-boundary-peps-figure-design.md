# Open-Boundary Toric-Code PEPS Figure Design

## Goal

Revise only the `plot_peps_graph` rendering so that every displayed toric-code
tensor has four virtual legs and four physical legs. The figure must show the
open boundary explicitly and preserve the isometric arrow convention in the
provided local-tensor reference. The stored `ToricCodePEPS`, its contracted
boundary vectors, normalization, graph topology, circuit model, and
wavefunction remain unchanged.

## Index and arrow convention

The directional order remains E/N/W/S. At each tensor,

```math
(v_E,v_N,v_W,v_S)=(\alpha,\beta,\gamma,\delta),
\qquad
(p_E,p_N,p_W,p_S)=(i,j,k,l).
```

The virtual isometry direction is

```math
\gamma \longrightarrow T \longrightarrow \alpha,
\qquad
\delta \longrightarrow T \longrightarrow \beta.
```

Consequently every horizontal virtual arrow points east and every vertical
virtual arrow points north. All four physical arrows point away from the
tensor along short northeast-sloping legs, matching the supplied local
diagram.

## Lattice rendering

Tensor `(row, col)` keeps the existing fixed position

```math
x=\mathrm{col}, \qquad y=R-\mathrm{row}+1.
```

The renderer constructs four virtual half-legs for every site. Two neighboring
half-legs meet to form an internal bond while remaining separately associated
with their tensors, so both local arrow directions are visible. A half-leg
without a neighbor ends as a dangling stub. Therefore a finite `R × C` figure
shows all `2R+2C` open-boundary virtual stubs even though the underlying PEPS
has already contracted those boundary indices with `|+\rangle` vectors.

Each tensor also owns four physical legs attached near its E/N/W/S virtual
ports. They are visually distinct from virtual legs and bonds but all slope
northeast, as in the reference diagram.

## Labels and visual hierarchy

The selected layout is the representative-site treatment:

- Every tensor displays all eight arrowed legs.
- Only the lower-left tensor `(R, 1)` carries the symbolic labels
  `α, β, γ, δ, i, j, k, l`.
- Other tensors retain identical arrow geometry without symbolic labels.
- Tensor nodes remain blue, virtual legs dark gray, and physical legs red.
- `show_index_labels=true` remains an optional diagnostic overlay for compact
  actual ITensor identifiers. It does not replace the representative symbolic
  labels.

This balances an explicit isometric convention with legibility on `3 × 3` and
larger publication figures.

## Public API and scope

The public signature is unchanged:

```julia
plot_peps_graph(peps::ToricCodePEPS; show_index_labels=false)
```

This is a rendering-only change. It must not modify `ToricCodePEPS`,
`peps_graph`, boundary contractions, physical-index ordering, Yao wire
allocation, circuit layers, postselection metadata, or normalization.

## Verification

Renderer tests will check:

- four virtual half-legs and four physical legs for every site;
- dangling virtual stubs on every open boundary;
- eastward horizontal and northward vertical virtual arrows;
- four northeast physical arrows at every tensor;
- symbolic labels only on the lower-left representative tensor;
- optional compact physical-index identifiers when requested;
- nonempty SVG and PDF vector export.

The complete existing package suite must remain unchanged and pass, providing
regression coverage for PEPS tensors, normalization, graph topology, Yao
circuits, postselection probabilities, and physical wavefunctions.
