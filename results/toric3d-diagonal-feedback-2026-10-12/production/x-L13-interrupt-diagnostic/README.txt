A deliberate Ctrl-C during X L13 p=0.04 batch_start=353 was misclassified as an invalid recovery by pre-fix code.
X replay records scalar triplets before incrementing completed_shots at the end of each 32-shot batch. The archived checkpoint therefore had completed_shots=352 while each joint-count row summed to 356.
Shots 353-356 were deterministically replayed from point seed d53a0848cbce9e52; their exact contributions were subtracted from the live checkpoint, whose rows now sum to completed_shots=352. The global and point invalid_recoveries counters were reset to zero, and the live failure TOML was removed after archival.
Joint-count contributions removed for logical_1/logical_2/state_failure were [1,1,1,1,0,0,0,0], [1,1,1,1,0,0,0,0], and [0,1,1,2,0,0,0,0].
Archived checkpoint SHA-256: c12ed2d845c5809aabe12765d5832b101329989fe434f4e2a3944840a2039a20
Archived failure SHA-256: 20a5f8b1d9f7092cfe45cdba09b652f9be89b2bb839c184430ec9ba9305c8404
Code fix: commit 1e4dbe9 rethrows InterruptException before recovery-failure bookkeeping.
The live sampling_seconds retains negligible sampling work from the interrupted batch; timing fields are therefore provenance totals rather than a pure completed-shot benchmark. Counts, RNG position, and scientific outcomes are unaffected.
