# Independent analysis implementation review

Status: **passed within the implementation-review scope**. The two identified
failure-handling issues are resolved. No remaining blocker was found in the
reviewed analysis wiring. Actual production numerical results have not yet been
reviewed; production was still running during this review.

Reviewed source bindings:

- `analyze.py`: `5542b748bcdd8db6e74cedc185f5902a10463723cfb0dc2a56a1a41800e3552f`
- `test_analysis.py`: `f93cfbd256522bca468993df3d80b77e8f6144eeae507ed59a6ad33e58be9c01`
- Frozen selection: `8fb02c7f5e476c0e9a5f69fb08209ef85088f7b576b23de07d18adc8b92a8750`

The review traced unequal-shot binomial likelihood weighting, the four-cell
Bonferroni anchor contrast, fresh-only X bootstrap inputs and frozen settings,
same-cell secondary pooling, all planned fit schedules, and the twelve Z profiles
with nuisance reoptimization. Old pooling independently resolves to 36 cells and
36,000 records: 12,000 X plus 24,000 Z. No fresh production outcomes were analyzed.

The failed free-parameter profile reproduction now returns an explicitly invalid
profile with a critical helper exception, empty confidence regions, and a note
that the unchanged helper's intermediate cache is unavailable. It no longer
aborts subsequent profiles. Numerical completion now checks missing/duplicate
results, optimizer and finite-value failures, profile disagreement/gaps, and
bootstrap failures. Boundary hits, nonmonotonicity and poor statistical fit remain
separate diagnostic flags.

Fresh verification reran seven focused failure/completion tests: **7 passed**
(exit code 0). An additional independent reproduction of the original failed
free-point path confirmed that no interval is emitted and the failure is retained.
The saved test log is `analysis-review-tests.log`, SHA-256
`a6f4aff2cc5e1bc6bc3c92af1f7222aa7c495aced60f0523a72562a51de205f2`.
The previously reported full 33-test result was not rerun for this focused review.

This review did not change analysis code, frozen numerical helpers, the selection,
simulation sources or production artifacts. Final artifact validation, independent
simulation replay and numerical result spot-checks remain subsequent work.
