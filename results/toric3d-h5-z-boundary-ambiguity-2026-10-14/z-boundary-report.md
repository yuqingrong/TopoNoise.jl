# H+5 Z-boundary ambiguity diagnostic

Status: `diagnostic_only`.

Logical Z is a state-independent diagnostic, not `|00_L>` infidelity. The decoder-visible information is the final star syndrome and observed upper parity.

| T | Lu | Lv | protocol | scope | locations | mixed groups | avoidable baseline failures | control / unavoidable floor |
|---:|---:|---:|:---|:---|---:|---:|---:|---:|
| 3 | 3 | 3 | `none` | final layer | 99 | 5 | 26 | 8 |
| 3 | 3 | 3 | `one_step` | final layer | 99 | 5 | 26 | 8 |
| 3 | 3 | 3 | `window_complete` | final layer | 99 | 5 | 26 | 8 |
| 3 | 3 | 3 | `full_prefix` | final layer | 99 | 5 | 26 | 8 |
| 5 | 5 | 5 | `none` | final layer | 275 | 9 | 117 | 14 |
| 5 | 5 | 5 | `one_step` | final layer | 275 | 9 | 117 | 14 |
| 5 | 5 | 5 | `window_complete` | final layer | 275 | 9 | 117 | 14 |
| 5 | 5 | 5 | `full_prefix` | final layer | 275 | 9 | 117 | 14 |
| 7 | 7 | 7 | `none` | final layer | 539 | 13 | 256 | 20 |
| 7 | 7 | 7 | `one_step` | final layer | 539 | 13 | 256 | 20 |
| 7 | 7 | 7 | `window_complete` | final layer | 539 | 13 | 256 | 20 |
| 7 | 7 | 7 | `full_prefix` | final layer | 539 | 13 | 256 | 20 |

The circuit-aware control removes the avoidable baseline mismatch but retains the information-theoretically unavoidable floor. It cannot establish a deterministic threshold for the current boundary instrument.

Counts assume equal single-location Z-fault priors; tied logical classes use deterministic order `00, 01, 10, 11`.
