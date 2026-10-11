# 3D toric-code preparation and feedback

Approved in conversation on 2026-10-11. Implement in the clean, dedicated `3dTC` branch.

## Scientific contract

Relative cubic slab, periodic y/z and rough x ends. Lx counts retained vertex planes; x edges have coordinates 0:Lx, transverse edges 1:Lx, and periodic coordinates 0:Ly-1, 0:Lz-1. Retain distinct cells at periodic length two. Include triangular boundary plaquettes. n=(3Lx+1)LyLz, k=1, dZ=Lx+1, dX=LyLz. Allow Lx>=1, Ly,Lz>=2.

Sequential vertices x/y/z increasing; outgoing x control receives H then CNOTs to incoming x, +y, -y, +z, -z. Ideal initialization, idle, measurement and correction. Independent X-only or Z-only Bernoulli p per active gate output. Compare final-only recovery with recovery after each x-layer on the entire prefix slab (moving triangular frontier). Future qubits cannot be corrected.

Record both logical Pauli diagnostics; |0_L> preparation failure is logical X only. Include logical-Z diagnostic lower bound (1-(1-2p)^(LyLz))/2 and its indistinguishable-fault proof, not a fidelity bound.

## Task 1: geometry, encoder, and fault engine

Own src/toric3d/geometry.jl, encoder.jl, noise.jl and test/toric3d/geometry.jl, encoder.jl, noise.jl. Implement ToricCode3D(Lx,Ly,Lz), star_checks, plaquette_checks, star_check_matrix, plaquette_check_matrix, existing generic qubit/logical support methods, toric3d_prefix(code,j) returning (code=prefix, qubits=local_to_full).

Toric3DEncoder contains code, operations::Vector{EncoderOperation}, checkpoints::Vector{Int}, prefixes::Vector{ToricCode3D}, prefix_qubits::Vector{Vector{Int}}. toric3d_encoder builds it. verify_toric3d_encoder checks ideal preparation and checkpoints algebraically.

Toric3DFaultRecord contains dimensions and x/z::BitMatrix (2 x number of operations, second H slot always false). sample_toric3d_faults(rng,encoder; p_x,p_z), with_toric3d_fault(encoder,record,step,qubit,pauli), propagate_toric3d_segment!(frame,encoder,record,first,last), toric3d_syndrome(code,frame) returning (stars,plaquettes), toric3d_logical_parities(code,frame) returning (logical_x,logical_z). Reuse PauliFrame and EncoderOperation, without changing 2D behavior.

## Task 2: decoders

Own src/toric3d/decoder.jl, test/toric3d/decoder.jl, CondaPkg.toml. Toric3DDecoderSettings keywords max_iter=50, bp_method="minimum_sum", ms_scaling_factor=1.0, schedule="parallel", omp_thread_count=1, osd_method="OSD_CS", osd_order=2. build_toric3d_decoder(code; p, x_backend=:bposd, settings=Toric3DDecoderSettings()). decode_toric3d(decoder,stars,plaquettes) returns PauliFrame recovery and requires syndrome consistency. toric3d_decoder_metadata returns backend names, versions, and settings (NamedTuple).

Z recovery uses PyMatching unit weights, full qubit output; X uses ldpc BP+OSD uniform p priors. :exact is minimum-weight data recovery for dimensions (1 or 2,2,2), enumerating the kernel of HZ (512 solutions for 2^3). This is not exact circuit maximum likelihood. Deterministic settings/ties, zero-syndrome fast path and p=0 safe. Fail invalid recovery. Tests compare exact weight/correctness, not arbitrary tied logical choices.

## Task 3: experiment and outputs

Root owns module integration, src/toric3d/scan.jl, plotting/output, test integration. compare_toric3d_feedback(rng; sizes=[3,5,7], error_rates=[0,.001,.003,.01,.02,.03,.05,.075,.1,.15], shots=1000, seed=1234, settings=..., x_backend=:bposd, output_dir=nothing) returns Toric3DComparison. sizes accepts integers for cubes or triples. Replay each sampled fault record in both protocols. Save completed points incrementally; CSV carries individual Wilson 95% intervals, paired differences and uncertainty, logical X/Z/any/state, both distances, full settings and versions. Save two-panel comparison PDF/SVG/PNG with analytic Z bound. No threshold fits.

## Task 4: documentation and runnable example

Own document/note.tex 3D section only, README appended 3D section, examples/compare_toric3d_feedback.jl. CLI calls compare_toric3d_feedback and save_toric3d_comparison(comparison,output_dir; basename="toric3d-feedback"). Flags --sizes 3,5,7; --error-rates comma-list; --shots 1000; --seed 1234; --output-dir path; --pilot-shots 10; --help. Time pilot on largest requested size at p=.03 before full scan, with independent pilot RNG. Default output results/toric3d-feedback-2026-10-11. Correct normalization, explicitly define geometry/schedule/feedback/noise, both metrics and decoder limitations, derive Z bound.

## Acceptance

Verify commutation/ranks/logicals and 28-qubit reference, prefixes/fresh controls, algebraic noiseless tableau, all eligible single faults, independent correction/frame equivalence, indistinguishable Z pair, decoder recovery syndrome, exact oracle, zero noise, reproducibility, and existing 2D regressions. Run timed pilot and full requested 60 paired scan points (120 protocol estimates), 1000 shots each. Review final changes and visually inspect plots.
