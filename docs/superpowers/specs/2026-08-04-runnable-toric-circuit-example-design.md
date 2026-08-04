# Runnable Toric-Code Circuit Example Design

## Goal

Provide a small command-line Julia example that lets a user construct an
executable toric-code PEPS circuit and export Yao's native circuit drawing
without creating a second Julia environment.

## Interface

Add `examples/generate_toric_circuit.jl`, invoked from the repository root:

```bash
julia --project=. examples/generate_toric_circuit.jl [rows] [cols]
```

`rows` and `cols` are optional positive integers and default to `2 2`. Invalid
or nonpositive dimensions produce a concise usage error and a nonzero exit.

## Behavior

The script will:

1. Construct `toric_code_peps(rows, cols)`.
2. Construct the reusable schedule with `sequential_circuit_graph(peps)`.
3. Construct the executable `Yao.ChainBlock` with `yao_circuit(circuit)`.
4. Print the lattice size, total qubit count, and diagonal site layers.
5. Create `examples/output/` relative to the script location.
6. Export the native Yao circuit to SVG, PDF, and PNG using a deterministic
   basename `toric_circuit_<rows>x<cols>`.
7. Print the absolute paths of the generated files.

The script constructs and draws the circuit but does not allocate or simulate
the exponentially large state vector. Simulation remains available to users
through the documented `zero_state` and `apply!` calls for suitably small
lattices.

## Repository Integration

- Add `/examples/output/` to `.gitignore` so generated figures are not staged.
- Add README instructions for `Pkg.instantiate()` and the example command.
- Use the root `Project.toml`; do not add a separate examples environment.
- Keep the existing Yao and CairoMakie plotting APIs unchanged.

## Validation

- Test the script with its default `2x2` dimensions in an isolated temporary
  output directory supplied through an internal `main` function argument.
- Verify the reported wire count is `4RC + R + C`, the printed layers begin at
  the lower-left site, and SVG/PDF/PNG files are nonempty.
- Test explicit `1x1` dimensions and invalid arguments.
- Run the complete package test suite and `git diff --check`.

