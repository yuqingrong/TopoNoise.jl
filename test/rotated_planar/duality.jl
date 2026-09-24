using Test, Random

# R is a clockwise quarter turn. Hadamard conjugation reverses each CNOT's
# arrows, not the temporal order: CNOT(c,t) -> CNOT(R(t),R(c)).
_duality_rotation(d) = [(column - 1) * d + d + 1 - row
                       for row in 1:d for column in 1:d]

@testset "As and Bp have matched full-lattice CNOT schedules" begin
    for d in (3, 5, 7, 9, 11, 13, 15), orientation in (:x_ns, :x_ew)
        code = RotatedPlanarCode(d; boundary_orientation=orientation)
        bp = rotated_planar_encoder(code; construction=:bp)
        as = rotated_planar_encoder(code; construction=:as)
        rotation = _duality_rotation(d)
        bp_blocks = filter(b -> b.kind === :plaquette, plaquette_blocks(bp))
        as_blocks = filter(b -> b.kind === :plaquette, plaquette_blocks(as))
        bp_cnots = [op.qubits for layer in gate_layers(bp)
                    for op in layer.operations if op.gate === :CNOT]
        as_cnots = [op.qubits for layer in gate_layers(as)
                    for op in layer.operations if op.gate === :CNOT]

        # These catch independently shelled families, a reversed time order,
        # or a local tree that only matches some plaquettes or distances.
        @test as_cnots == [[rotation[t], rotation[c]] for (c, t) in bp_cnots]
        @test [sort(b.source_support) for b in as_blocks] ==
              [sort(rotation[b.source_support]) for b in bp_blocks]
        @test getproperty.(as_blocks, :representative) ==
              rotation[getproperty.(bp_blocks, :representative)]
        @test Set(Tuple(sort(rotation[s])) for s in a_s_checks(code)) ==
              Set(Tuple(sort(s)) for s in b_p_checks(code))
        @test Set(Tuple(sort(rotation[s])) for s in b_p_checks(code)) ==
              Set(Tuple(sort(s)) for s in a_s_checks(code))
        @test as.logical_state === :plus
        @test bp.logical_state === :zero
        @test verify_encoder_tableau(as)
        @test verify_encoder_tableau(bp)

        # Native As still uses direct product-input preparation, with no
        # logical-parity seed or basis-conversion circuit.
        as_h = [only(op.qubits) for layer in gate_layers(as)
                for op in layer.operations if op.gate === :H]
        @test sort(as_h) == sort(rotation[getproperty.(bp_blocks, :representative)])
        @test all(b.source_check ∉ (:logical_parity, :logical_conversion)
                  for b in plaquette_blocks(as))
    end
end

function _duality_failures(code, decoders, frames)
    syndromes = measure_syndrome.(Ref(code), frames)
    a_s = BitMatrix([s.a_s[j] for s in syndromes, j in eachindex(first(syndromes).a_s)])
    b_p = BitMatrix([s.b_p[j] for s in syndromes, j in eachindex(first(syndromes).b_p)])
    prediction = decode_logical_parities(decoders, a_s, b_p)
    actual_x = [isodd(count(f.x[logical_z_support(code)])) for f in frames]
    actual_z = [isodd(count(f.z[logical_x_support(code)])) for f in frames]
    return (x=xor.(prediction.logical_x, actual_x),
            z=xor.(prediction.logical_z, actual_z))
end

@testset "Matched H, CNOT, and plaquette faults exchange decoded X and Z" begin
    for d in (3, 9, 11, 13, 15), orientation in (:x_ns, :x_ew)
        code = RotatedPlanarCode(d; boundary_orientation=orientation)
        rotation = _duality_rotation(d)
        bp = rotated_planar_encoder(code; construction=:bp)
        as = rotated_planar_encoder(code; construction=:as)
        bp_clean = sample_fault_record(MersenneTwister(1), bp, CircuitPauliNoise(0))
        as_clean = sample_fault_record(MersenneTwister(1), as, CircuitPauliNoise(0))
        bp_frames, as_frames = PauliFrame[], PauliFrame[]

        # Exhaust X/Z faults on every native H and both CNOT operands.
        for (bp_index, as_index) in zip(eachindex(bp_clean.steps), eachindex(as_clean.steps)),
            q in bp_clean.steps[bp_index].support, (pauli, dual) in ((:X, :Z), (:Z, :X))
            push!(bp_frames, propagate_pauli_frame(bp,
                with_pauli_fault(bp_clean, bp_index, q, pauli)))
            push!(as_frames, propagate_pauli_frame(as,
                with_pauli_fault(as_clean, as_index, rotation[q], dual)))
        end

        # Couple complete many-fault histories, including noisy H gates,
        # for both clocks. Independent RNG streams need not match shotwise.
        rng = MersenneTwister(932 + d)
        for clock in (:gate_layer, :plaquette), _ in 1:32
            as_clean = sample_fault_record(rng, as, CircuitPauliNoise(0; clock))
            bp_faults = sample_fault_record(rng, bp,
                CircuitPauliNoise(0; p_x=0.13, p_z=0.21, clock))
            steps = CircuitFaultStep[]
            for (source, target) in zip(bp_faults.steps, as_clean.steps)
                bits = Dict(rotation[q] => (z, x)
                            for (q, x, z) in zip(source.support, source.x, source.z))
                push!(steps, CircuitFaultStep(target.step_index, target.block_kind,
                    target.block_order, target.support,
                    [bits[q][1] for q in target.support],
                    [bits[q][2] for q in target.support]))
            end
            as_faults = CircuitFaultRecord(clock, d, orientation, :as, :plus, steps)
            push!(bp_frames, propagate_pauli_frame(bp, bp_faults))
            push!(as_frames, propagate_pauli_frame(as, as_faults))
        end

        @test all(a.x[rotation] == b.z && a.z[rotation] == b.x
                  for (a, b) in zip(as_frames, bp_frames))
        decoders = build_matching_decoders(code)
        bp_failures = _duality_failures(code, decoders, bp_frames)
        as_failures = _duality_failures(code, decoders, as_frames)
        @test as_failures.z == bp_failures.x
        @test as_failures.x == bp_failures.z
    end
end
