# TN experiment with the uploaded figure's preparation

This experiment uses the integrated `TopoNoise` boundary-MPS decoder. It adds
the complete initial seed repair present in the archived no-feedback figure;
the production decoder in `src/` is unchanged.

## Physical experiment and decoder information

1. Apply independent Z faults with probability `p` to every seed-H output.
2. Measure the seed checks and apply the existing unit-weight seed MWPM repair.
3. Run the diagonal H+6CNOT preparation with independent Z faults on every
   gate-output qubit. Apply no further preparation feedback.
4. Measure the final syndrome and infer the two-bit logical error class.

The TN decoder retains both the seed and final syndrome records and the known
seed recovery. The paired final-MWPM baseline uses the final syndrome only.
The physical preparation is matched to the archived figure; the information
used for the final decision differs. All measurements and corrections here
are ideal. Failure means either logical bit is inferred incorrectly.

## Why the repaired-seed likelihood factors

Let `a` and `b` be independent raw seed and transfer faults, and `R(s0)` the
actual seed recovery. The repaired seed has zero final syndrome, so the final
syndrome constrains only `b`. In the common logical-mask order `(00,10,01,11)`,

```text
q(ell | s0, s) = sum_a q_seed(a | s0)
                         q_transfer(ell XOR a XOR B_seed*R(s0) | s).
```

Here the sum is over four *logical masks*, after summing the underlying raw
faults in each TN. `model.jl` builds the two topologies by splitting raw columns
before combining equivalent faults. It retains signed numerical weights and
requires both component contractions to be valid. There is no exact fallback.

## Run and resume

Use the same Julia project and Python environment as `TopoNoise`, one Julia
thread and one BLAS thread. Arguments are:

```text
julia --project=. --startup-file=no --threads=1 examples/tnml_repaired_seed/run.jl \
    L p shots chi rng_seed output_dir exact=true
```

`exact=true` is optional and supported only for `L <= 3`. It independently
contracts a dense, untruncated detector frontier for each component. The exact
comparison does not replace the TN decision. Cutoff is `1e-18`; the numerical
negative-weight tolerance is `1e-10`.

The saved four-point pilot can be resumed sequentially with
`python examples/tnml_repaired_seed/run_pilot.py --root results/tnml-matched-figure-2026-10-27`.
Use `--julia /absolute/path/to/julia` if Julia is not on the executable search
path. Its frozen `pilot-selection.json` specifies 20 independent L3 histories
at each of four rates, with an exact reference for every shot. Run
`plot.py --root results/tnml-matched-figure-2026-10-27` to regenerate the figures
and machine-readable counts. Diagnostic histories are excluded by default.

The runner saves every requested raw history before attempting contractions.
It records source hashes, runtime versions, syndromes, recovery, truth used
only for scoring, TN results, exact comparisons, and elapsed times. Resuming
the identical command skips saved outputs and rejects changed source,
configuration, or histories. Invalid/exception outputs remain recorded;
use a fresh output directory to evaluate different settings on the same RNG
seed.

`supervise.py` optionally wraps the command with per-contraction and total
time limits plus a resident-memory limit. It stops only its own worker process
group. When process inspection is unavailable, it reports that the memory
limit was not enforced; time limits remain active.

```text
python examples/tnml_repaired_seed/supervise.py \
    --output output_dir --decode-seconds 180 --total-seconds 900 --rss-gib 4 \
    -- julia --project=. --startup-file=no --threads=1 \
       examples/tnml_repaired_seed/run.jl L p shots chi rng_seed output_dir
```

## Interpretation and checks

Finite bond dimension gives an approximate ML decoder. Numerical validity and
small local SVD truncation errors do not certify posterior convergence.
Compare exact L3 logical probabilities and decisions; check larger bond
dimensions before making larger-size threshold claims.

For `N` requested histories, `F` known TN failures and `U` unknown outcomes,
the realized sample failure fraction is only identified within
`[F/N, (F+U)/N]`. These are missing-outcome bounds, not confidence intervals.
Never divide by just the valid shots. Complete binary-outcome points can use
95% Wilson intervals; small pilot samples do not establish a threshold.

Run `test_model.jl` for physical replay, exhaustive L2 seed, and independent
probability-oracle checks. Run `test_runner.jl` for checkpoint, identity,
truth separation, tie-aware comparison, and missing-outcome accounting checks.
