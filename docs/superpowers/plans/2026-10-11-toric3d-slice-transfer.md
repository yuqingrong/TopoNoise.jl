# 3D slice-transfer and matching feedback implementation plan

> **For agentic workers:** Use subagent-driven-development for disjoint components, with root integration and independent scientific review.

**Goal:** Implement and test the user's approved slice-transfer construction and active-slice feedback while retaining independent Bernoulli X-only/Z-only faults after every H/CNOT output.

**Architecture:** Add a separate diagonal-lattice state-preparation subsystem. An old slice controls a fresh GHZ triad at each vertex; a measured six-edge compatibility check defines a triangular MWPM graph. Matching routes determine ideal Pauli corrections on new targets only. Preserve the existing rough-slab implementation and results.

**Tech stack:** Julia, binary Pauli/stabilizer algebra, PythonCall/PyMatching, CairoMakie, TOML.

**Spec:** User's approval to try steps 1 and 2 in the preceding paper comparison, explicitly preserving gate-output noise. This is a derived 3D realization of arXiv:2609.13383, not an author-supplied circuit or threshold reproduction. The interface decisions below concretize the approved work.

## Global constraints

- Ideal initialization, measurements, feed-forward Paulis, and idles. Faults occur independently after each participating H/CNOT output, separately for X-only and Z-only channels.
- Fixed H/CNOT schedule and pre-sampled faults, shared between unassisted and feedback runs.
- No decoder receives the sampled faults. Every full matching recovery must reproduce its syndrome; routed partial corrections report their residual syndrome explicitly.
- Initial slice is all-zero; the resulting boundary stabilizers specify a unique state (k=0). Do not reuse old rough-slab logical operators/distances or interpret zero final syndrome as a threshold test.
- Default routing moves each paired endpoint one edge toward its partner per slice. Expose full-path routing as a separately labeled stronger ideal-feedback option.
- No commits, staging, or destructive changes; existing uncommitted work and generated note files belong to the user/current task.

## Mathematical contract

Coordinates use t=x+y+z, u=x, v=y, with periodic u,v. Incoming edge E[t,d,u,v] ends at vertex(t,u,v), d=1,2,3; spatial shifts s1=(1,0),s2=(0,1),s3=(0,0). Slices run 0:steps. For new vertex r,

Kij(b,r)=b_j(r-si) xor b_i(r-sj), Q(b,r)=K12 xor K13 xor K23.

The gadget is H(new1), CNOT(new1,new2), CNOT(new1,new3), the two old K12 controls into new2, then the two old K13 controls into new3. This has seven gates and thirteen active outputs. A measured old-slice recovery C changes only new targets by F(C)=(0,K12(C),K13(C)). With no new gate faults, Qnew(r)=[Qold xor Q(C)](r-s1).

Q has N=LyLz rows and 3N columns, column weight two and rank N-1. Keep all physical edge IDs, including periodic length two. Physical plaquettes have new edges i,j and their two Kij controls. Bulk stars have three incoming plus three outgoing edges, and the final frontier has three-edge X stars. Initial boundary checks are individual Z on E0. Ranks: n=3N(steps+1), rankX=N steps, rank(plaquettes)=(2steps+1)N-1, rank(fullZ)=(2steps+3)N. Every single Z gate-output fault has final weight at most three.

## Review focus

1. Thin/rectangular/periodic-length-two lattices: distinct physical edges and correct ranks.
2. Intermediate measurements: Q must be an actual prepared-state stabilizer, not an assumed future plaquette.
3. Partial routing: clear separation between full valid recovery and one-step correction with remaining defects.
4. New faults during a transfer: do not assert outgoing syndrome is always cleared by a pre-transfer decision.
5. Physics interpretation: global classical matching and parallel Pauli gates are not a depth-one local channel; report propagation/cluster diagnostics only, with no threshold claim.

## Tasks

### 1. Diagonal lattice and transfer circuit

Files: src/toric3d/slice_geometry.jl, slice_encoder.jl; test/toric3d/slice_geometry.jl, slice_encoder.jl.

Interfaces: Toric3DSliceLattice(steps,Ly,Lz), fields steps,Ly,Lz. Generic data_qubit_count; data_qubit_index(code,d::Integer,t,u,v) validates explicit coordinates, direction1:3; data_qubit_coordinate returns (direction,t,u,v). toric3d_slice_qubits(code,t) returns full indices in direction,u,v order. Generic star_checks, plaquette_checks and corresponding matrices; toric3d_slice_boundary_checks(code); toric3d_slice_checks(code,t); toric3d_slice_check_matrix(code) is local N×3N. _toric3d_slice_pair(code,t,u,v,i,j) returns two global old indices used by Kij for new time t.

Toric3DSliceEncoder fields code,operations,checkpoints; toric3d_slice_encoder(code); verify_toric3d_slice_encoder(encoder) tests all prefixes with signed binary conjugation.

- [x] Write and run failing API/algebra tests.
- [x] Implement geometry and circuit; verify commutation/ranks, actual H/CNOT order, initial boundary and every checkpoint Q/star/plaquette eigenvalue.
- [x] Report targeted verification; root integrates module includes.

### 2. Matching, feed-forward and faults

Files: src/toric3d/slice_feedback.jl; test/toric3d/slice_feedback.jl. Root owns.

Interfaces: build_toric3d_slice_decoder(code), decode_toric3d_slice(decoder,syndrome;routing=:one_step) returns full recovery, applied route, pairs and residual syndrome. Extend sample_toric3d_faults/with_toric3d_fault/propagate_toric3d_segment! for Toric3DSliceEncoder using existing Toric3DFaultRecord (Lx stores steps). replay_toric3d_slice(encoder,faults,decoder;feedback=true,routing=:one_step,record_history=false) returns frame, measured frontier histories and physical syndromes.

- [x] Write failing tests for even syndrome, recovery consistency, shortest paths, odd-syndrome rejection, one-step endpoint motion and length-one paths.
- [x] Implement unit-weight matching and deterministic shortest paths; validate full recovery.
- [x] Implement target-only feed-forward and fixed-record replay.
- [x] Verify exhaustive single faults (especially Z support≤3), virtual-control equivalence, noiseless behavior, seeded replay, and no sampled-fault access by decoder.

### 3. Pilot and delivery

Files: src/toric3d/slice_scan.jl; examples/compare_toric3d_slice_feedback.jl; README addition; results/toric3d-slice-feedback-2026-10-11.

- [x] Configurable paired scan, default pilot sizes3,5,7, six p values0,.001,.003,.01,.03,.05,200 shots/channel,seed1234; preserve original full-scan files.
- [x] Record final bulk star/plaquette syndrome densities, final compatibility syndrome fraction, and longitudinal extent of largest connected plaquette-syndrome cluster normalized by depth. Save paired aggregate statistics, settings and noise/gate metadata incrementally, plus CSV/PDF/SVG/PNG.
- [x] Explain graph/circuit derivation, k=0 boundaries, one-step versus full-path routing, and why diagnostics do not establish a threshold.
- [x] Independent review, targeted tests plus existing regression suite, output checks and visual inspection. Full `test/runtests.jl` exited successfully; the new slice subsystem passed 2,971 assertions. The 36-point scan completed, and its statistics and rendered figure were independently checked.
