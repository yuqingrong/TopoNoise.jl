```@meta
CurrentModule = TopoNoise
```

# Isometric planar-code capacity

[`IsometricPlanarCodeModel`](@ref) defines an open planar logical-code
experiment parameterized by its direct distance `d`. It prepares ``|0_L\rangle``
with west/east ``|+\rangle`` and south/north ``|0\rangle`` virtual boundaries.
Independent `X` errors occur only on internal virtual carriers; preparation,
local gates, and parity checks are ideal.

The decoder receives the ``(d-1)^2`` plaquette checks and one local north and
south check for every horizontal boundary bond. The held-out logical frame is
the parity of residual west-boundary vertical bonds, making a west-to-east
virtual string the shortest odd logical path, of weight `d`.

```julia
using TopoNoise, Random

model = IsometricPlanarCodeModel(5)
point = estimate_isometric_planar_capacity(
    MersenneTwister(1234), model, 0.10; shots=10_000, batches=100)
```

[`scan_isometric_planar_capacity`](@ref) evaluates a distance and error-rate
grid; [`estimate_isometric_planar_crossings`](@ref) reports adjacent-distance
crossings. `d=2` is a supported diagnostic geometry with a fixed tie rule, but
is excluded from finite-size threshold fitting. For a production scaling scan,
use at least three non-diagnostic distances, for example:

```bash
julia --project=. examples/scan_isometric_planar_threshold.jl \
  --distances 7,9,11,13,15 \
  --output-dir results/isometric-planar-code/d7-15
```

The script writes raw scan and crossing CSVs, a capacity figure, and scaling
diagnostics. [`fit_isometric_planar_scaling`](@ref) fits a cubic master curve
to logical-failure rates in the configured fit window under
``P_{\mathrm{fail}} = F((p-p_c)d^{1/\nu})``. Use
[`diagnose_isometric_planar_scaling`](@ref) to inspect the fit, bootstrap
refits, loss surface, optimizer path, batch variation, and leave-one-distance
or fit-window sensitivity checks.

[`isometric_planar_encoder`](@ref) keeps its terminal carriers as physical
outputs, so it has ``4d^2`` output wires. The small-distance
[`sample_isometric_planar_yao_trajectory`](@ref) path is a detector reference;
the estimator computes detector parities algebraically and decodes with
PyMatching without allocating an exponential state vector.
