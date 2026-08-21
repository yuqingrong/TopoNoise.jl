# Generated results

Runnable generators live in `examples/`; their default artifacts are saved in
the experiment-specific folders below. Run these commands from the repository
root.

| Result folder | Generator |
| --- | --- |
| `results/toric-circuit/` | `julia --project=. examples/generate_toric_circuit.jl 2 2` |
| `results/toric-trajectories/` | `julia --project=. examples/scan_toric_trajectories.jl --p-min 0.35 --p-max 0.65 --p-step 0.01` |
| `results/rotated-code-capacity/` | `julia --project=. examples/scan_rotated_code_capacity.jl --p-min 0.05 --p-max 0.15 --p-step 0.01` |
| `results/isometric-planar-code/` | `julia --project=. examples/scan_isometric_planar_threshold.jl` |

The scan examples accept `--output-dir PATH` for an alternate destination.
CSV data can be versioned; rendered PDF, PNG, and SVG figures are ignored.
