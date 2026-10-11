# 3D toric preparation: final-only decoding versus layer feedback

Completed on 2026-10-11 with seed 1234: cubic sizes 3, 5, 7; physical
probabilities 0, 0.001, 0.003, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15;
1000 paired shots for each size, probability, and noise channel. All 60
channel points completed, comprising 60,000 fault records replayed under
both protocols (120,000 protocol realizations). Every recovery reproduced
its input syndrome. No threshold was fitted.

![Paired logical-error comparison](toric3d-feedback.png)

## Files

- [CSV](toric3d-feedback-raw.csv): 240 rows, one per point and metric,
  containing both protocol counts, rates, Wilson intervals, paired
  discordant counts, differences, uncertainties, decoder settings, and versions.
  These are aggregate statistics, not individual fault histories.
- [PDF](toric3d-feedback.pdf), [SVG](toric3d-feedback.svg), and
  [PNG](toric3d-feedback.png): the same two-panel comparison, including the
  analytic lower bound for the logical-Z diagnostic.
- [Metadata](toric3d-feedback-metadata.toml): complete experiment settings
  and provenance; `complete = true` and 60 completed channel points.
- [Run log](run.log): warmup, timed pilot, and completed-point timings.

## Selected results

At **L = 7, p = 0.10**, X-only noise gave:

| Protocol | Logical-X / preparation failure | 95% Wilson interval |
| --- | ---: | ---: |
| Final-only decoding | 50.8% | 47.7–53.9% |
| Layer feedback | 12.1% | 10.2–14.3% |

The paired feedback-minus-final difference is **−38.7 percentage points**,
with an approximate conservative 95% interval of **−43.6 to −33.4 points**.
Feedback also reduced logical-X failure at several other scanned points;
the complete CSV retains all points, including those without a resolved
difference. Zero observed failures at 1000 shots still have a Wilson upper
limit of approximately 0.383%.

For Z-only noise at the same size and probability, the logical-Z diagnostic
was 51.0% for final-only and 48.4% for feedback. The paired interval for
their difference includes zero. The theoretical diagnostic lower bound is
approximately 0.499991. A finite-sample estimate can fall below this bound;
the bound constrains the underlying probability, not every Monte Carlo
estimate. Logical Z leaves the prepared logical-zero state unchanged, so
these Z diagnostic rates are not preparation infidelities. Both residual
Pauli rates, their union, and state failure are recorded for every point.

## Model and statistical conventions

The slab has retained vertex planes `x=1:Lx`, dangling x edges at both rough
ends, and periodic y/z directions. Each gadget applies H to a fresh outgoing
x edge, then CNOTs to the incoming x, +y, −y, +z, −z edges. Every active gate
output independently receives a rate-p fault in the selected X-only or
Z-only channel. Initialization, idles, measurement, and correction are ideal.
Layer feedback corrects the entire prepared prefix, including its truncated
frontier, before subsequent gates. The two protocols use identical faults.

Z recovery uses unit-weight PyMatching 2.4.0. X recovery uses ldpc 2.4.1
BP+OSD with uniform physical-p priors, minimum-sum updates, 50 iterations,
scaling factor 1, parallel updates, one thread, and OSD_CS order 2. These
data-error priors neglect correlations induced by the preparation circuit.
The code provides a separate exact minimum-data-weight small reference;
this production study uses BP+OSD, not exact circuit-level maximum likelihood.

Individual intervals are 95% Wilson score intervals. Paired differences use
feedback-only and final-only discordant counts; their standard errors use
the paired sample variance. The reported difference interval combines two
97.5% Wilson discordance intervals with a Bonferroni construction, giving an
approximate conservative 95% interval, including when no discordance occurs.

## Reproduction and validation

Run from the repository root:

```bash
julia --project=. examples/compare_toric3d_feedback.jl
```

The example's default arguments reproduce this grid. The Julia API also
accepts rectangular dimensions, arbitrary scan probabilities, shot counts,
an explicit RNG, and decoder controls; see the [project README](../../README.md)
and [runnable example](../../examples/compare_toric3d_feedback.jl).

This machine used Julia 1.12.6 and the existing managed Python 3.12.13.
CondaPkg's mapping service failed during environment resolution, so this
successful run selected the managed interpreter explicitly:

```bash
JULIA_CONDAPKG_BACKEND=Null \
  JULIA_PYTHONCALL_EXE="$PWD/.CondaPkg/.pixi/envs/default/bin/python" \
  julia --project=. examples/compare_toric3d_feedback.jl
```

Warmup, including compilation, took 11.92 seconds. The independent pilot
used 10 paired shots per channel at L=7, p=0.03 and took 1.03 seconds.
Its largest-size extrapolation was conservative: summed sampling and
decoding times for all completed points were 495.83 seconds (8.26 minutes),
excluding decoder construction, startup, and plotting. Completed points
were saved incrementally.

The full existing 2D and new 3D regression suite passed. Checks include
stabilizer ranks and commutation, periodic length two, logical operators,
noiseless prefix preparation, fresh controls, physical/frame correction
equivalence, every small-code gate-output single fault, indistinguishable
Z-fault pairs, exact-reference consistency, BP+OSD validity, and seeded
reproducibility. An independent audit recomputed the saved statistics and
confirmed the requested grid, settings, zero-noise success, and logical
metric identities. The updated [3D note](../../document/note.tex) compiles;
the figure exports were rendered and visually checked.

The note is a local file already excluded by the repository's `*.tex`
ignore rule; preserve it explicitly when transferring these changes.
