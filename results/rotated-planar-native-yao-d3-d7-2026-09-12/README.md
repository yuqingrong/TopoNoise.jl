# Native Yao encoder drawings

These diagrams were exported from the current `yao_encoder(encoder)` using
`Yao.vizcircuit`, not reconstructed from the reference plaquette drawings.

- Constructions: `:as` (Z-type checks) and `:bp` (X-type checks).
- Distances: 3, 5, 7; respectively 9, 25, 49 data-qubit wires.
- Prepared state: `logical_state=:zero`, with `boundary_orientation=:x_ns`.
- Ideal encoders only: no sampled errors or final measurement circuit is drawn.
- Each encoder is kept as a PDF: `as-d3-zero-x_ns.pdf`, etc.
- The PDFs preserve the white backgrounds added to the native SVG canvases
  for dark-mode readability. No gate symbols, positions, or connections
  were changed during conversion.

## Reading the diagrams

`qN [row,column] |0>` identifies a data qubit, its row-major lattice coordinate,
and its initial physical state. All wires start in physical |0>; the displayed
H gates supply the required superposition preparation. H denotes Hadamard.
A filled dot is a CNOT control and a circled plus is its target. A vertical
line merely passing across another wire does not couple to that wire.

Yao uses a compact layout: independent gates can share a displayed column.
The drawing preserves the quantum circuit, but horizontal position is not the
literal global gate-by-gate noise clock. The stored simulation schedule was
not changed for drawing. The As logical-sector preparation includes its actual
parity-preparation CNOTs; no SWAP routing or additional gates were inserted.

## Checked encoders

| d | Construction | Data qubits | H | CNOT | Total gates |
|---|---|---|---|---|---|
| 3 | As | 9 | 4 | 10 | 14 |
| 3 | Bp | 9 | 4 | 8 | 12 |
| 5 | As | 25 | 12 | 31 | 43 |
| 5 | Bp | 25 | 12 | 28 | 40 |
| 7 | As | 49 | 24 | 63 | 87 |
| 7 | Bp | 49 | 24 | 60 | 84 |

All six schedules passed `verify_encoder_tableau`: every As and Bp stabilizer
has eigenvalue +1, and logical Z has eigenvalue +1. For d=3, both constructions
were also checked by exact Yao statevector execution and observable expectations.
The Yao blocks were checked against the complete stored elementary operation
sequence. There were 36 passing encoder/observable checks and 24 passing
initial export checks. Plotting did not allocate large-d statevectors.

## Draw another size directly

From a Julia session using this repository's project environment:

```julia
using TopoNoise
import Yao

d = 9                       # choose an odd integer >= 3
construction = :as          # or :bp
code = RotatedPlanarCode(d; boundary_orientation=:x_ns)
encoder = rotated_planar_encoder(code; construction, logical_state=:zero)
@assert verify_encoder_tableau(encoder)
Yao.vizcircuit(yao_encoder(encoder);
    format=:pdf, filename="$(construction)-d$(d).pdf")
```

The PDF output preserves vector graphics for zooming larger circuits. Changing `d`
changes the actual encoded lattice, not merely the drawing scale.

Yao 0.9.3; YaoPlots 0.9.6. Drawing API documentation:
https://docs.yaoquantum.org/latest/man/plot.html

Encoder source SHA-256 at export:
`b74156e7a3510b60341d867e68ede001cbf850d238327035a02c1a28d12b43e3`.
Geometry source SHA-256 at export:
`6f82c3bd05da708489ba6bb70912093a7b5be2e0765588f5cb133beb02521749`.
