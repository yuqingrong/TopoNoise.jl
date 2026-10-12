A deliberate Ctrl-C stopped the original Z L21 ten-rate shard after p=0.04 completed and p=0.045 had completed 32 shots. The process had been launched before commit 1e4dbe9 loaded the InterruptException fix, so it misclassified the stop at p=0.045 batch_start=33 as an invalid recovery.

The archived checkpoint records completed_shots=[2000,2000,2000,2000,2000,2000,32,0,0,0]. Every joint-count row sums to its completed_shots value, so the interrupted batch contributed no scored triplets. The six p<=0.04 points are complete and clean. The partial main-shard p=0.045 point is excluded from analysis; an independently derived exact-seed shard supplies that point.

The live checkpoint was repaired only by resetting the global and p=0.045 invalid_recoveries counters to zero. The failure TOML was removed from the live resumable directory after archival. Timing fields retain work from the interrupted batch and are provenance totals rather than pure completed-shot benchmarks. Counts, RNG position, and scientific outcomes are unchanged.

Archived checkpoint SHA-256: e5f499681b96490a60b8325c29ca98d8a5427d82c350c001d9511ae757b40a1a
Archived failure SHA-256: 6df348f36d9b0608fa323014fef5cc5fc98b9e944863e9d9aaff60b2a46c372f
Repaired live checkpoint SHA-256: 0f2f90d5887a826f65de5a89b9091c02fe75f623659fec3681f29a37859b5fc7
Code fix: commit 1e4dbe9 rethrows InterruptException before recovery-failure bookkeeping.
