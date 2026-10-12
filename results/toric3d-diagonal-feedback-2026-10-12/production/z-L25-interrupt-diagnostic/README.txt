A deliberate Ctrl-C during Z L25 p=0.03 batch_start=385 was misclassified as an invalid recovery by pre-fix code.
The archived checkpoint records completed_shots=384 and all three joint-count rows sum to 384; the interrupted batch contributed no scored triplets.
The live checkpoint was repaired only by resetting the global and point invalid_recoveries counters to zero. The failure TOML was removed from the live resumable directory after archival.
Archived checkpoint SHA-256: 5eedc686f7e8d0fbae6ba62c4db1f7c163de11e0c49683fd371a7bdaa7e6845b
Archived failure SHA-256: 61a24b64ea021f042cefcc26d21836b1cc67d9aa182c8fecc9ea60a4f5983b68
Code fix: commit 1e4dbe9 rethrows InterruptException before recovery-failure bookkeeping.
The live sampling_seconds and replay_seconds retain work from the interrupted batch; timing fields are therefore provenance totals rather than a pure completed-shot benchmark. Counts, RNG position, and scientific outcomes are unaffected.
