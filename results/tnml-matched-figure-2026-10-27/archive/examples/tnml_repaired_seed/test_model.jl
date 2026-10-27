using Test, Random, SparseArrays, TopoNoise

if isfile(joinpath(@__DIR__, "model.jl"))
    include("model.jl")
else
    @eval module RepairedSeedTN end
end

using .RepairedSeedTN

# Independent oracles use physical supports and ordinary GF(2) probability
# updates, never reduced TN edges or combine_posteriors.
bits_of(mask, n) = BitVector(!iszero(mask & (1 << (k-1))) for k in 1:n)
bits_mask(bits) = sum(Int(bit) << (k-1) for (k, bit) in enumerate(bits); init=0)
parities(checks, frame) = BitVector(isodd(count(frame[check])) for check in checks)

function raw_signatures(ctx, columns)
    H = star_checks(ctx.code)
    X = toric3d_encoded_logical_x(ctx.code)
    n = length(H)
    return [begin
        physical = Bool.(Vector(ctx.model.matrix[:, column]))
        bits_mask(parities(H, physical)) | (bits_mask(parities(X, physical)) << n)
    end for column in columns]
end

function transfer_dp(initial, signatures, p)
    distribution = copy(initial)
    scratch = similar(distribution)
    for signature in signatures
        @inbounds for state in 0:length(distribution)-1
            scratch[state+1] = (1-p)*distribution[state+1] +
                p*distribution[xor(state, signature)+1]
        end
        distribution, scratch = scratch, distribution
    end
    return distribution
end

function physical_record(ctx, raw)
    seed = raw[ctx.seed_columns]
    frame = BitVector(isodd.(ctx.model.matrix * UInt8.(raw)))
    seed_frame = falses(data_qubit_count(ctx.code))
    seed_frame[ctx.seed_qubits] .= seed
    seed_syndrome = parities(toric3d_encoded_seed_checks(ctx.code).x, seed_frame)
    seed_recovery = decode_toric3d_z_feedback(ctx.seed_decoder, seed_syndrome)
    repair_frame = falses(length(frame))
    repair_frame[ctx.seed_qubits] .= seed_recovery
    frame .⊻= repair_frame
    recovery_mask = bits_mask(parities(toric3d_encoded_logical_x(ctx.code), repair_frame))
    final_syndrome = parities(star_checks(ctx.code), frame)
    true_mask = bits_mask(parities(toric3d_encoded_logical_x(ctx.code), frame))
    return (;seed_syndrome, seed_recovery, recovery_mask, final_syndrome, true_mask, frame)
end

# Catches missing recovery shift, swapped logical-bit order, and clipping.
@testset "four-sector XOR labels and signed weights" begin
    for a in 0:3, b in 0:3, recovery in 0:3
        seed, transfer = zeros(4), zeros(4)
        seed[a+1] = transfer[b+1] = 1
        expected = zeros(4)
        expected[xor(xor(a, b), recovery)+1] = 1
        @test combine_posteriors(seed, transfer, recovery) == expected
    end
    signed = [-1e-12, 1+1e-12, 0., 0.]
    @test combine_posteriors(signed, [1.,0,0,0], 0) == signed
    @test_throws DimensionMismatch combine_posteriors(ones(3), zeros(4), 0)
    @test_throws ArgumentError combine_posteriors([NaN,0,0,1], [1.,0,0,0], 0)
    @test_throws ArgumentError combine_posteriors([1.,0,0,0], [1.,0,0,0], 4)
end

@testset "physical seed repair and independent exact posterior" begin
    ctx = build_context(2)
    H, X = star_checks(ctx.code), toric3d_encoded_logical_x(ctx.code)
    H0 = toric3d_encoded_seed_checks(ctx.code).x
    n, n0 = length(H), length(H0)
    nseed = length(ctx.seed_columns)
    @test nseed == 12
    @test length(ctx.transfer_columns) == 104
    @test n == 16
    @test ctx.seed_decoder.window.qubits == ctx.seed_qubits
    @test sort(vcat(ctx.seed_columns, ctx.transfer_columns)) == collect(axes(ctx.model.matrix, 2))

    # Enumerate all 4096 seed histories before adding any transfer faults.
    histories = NamedTuple[]
    recovery_by_syndrome = Dict{Int,BitVector}()
    nonzero_recovery_seed = nothing
    loop_seeds = Dict{Int,Int}()
    for rawmask in 0:(1 << nseed)-1
        seed = bits_of(rawmask, nseed)
        physical = falses(data_qubit_count(ctx.code))
        physical[ctx.seed_qubits] .= seed
        s0 = parities(H0, physical)
        s0mask = bits_mask(s0)
        recovery = get!(recovery_by_syndrome, s0mask) do
            decode_toric3d_z_feedback(ctx.seed_decoder, s0)
        end
        repair = falses(length(physical))
        repair[ctx.seed_qubits] .= recovery
        rawlogical = bits_mask(parities(X, physical))
        recoverylogical = bits_mask(parities(X, repair))
        physical .⊻= repair
        @test !any(parities(H0, physical))
        @test !any(parities(H, physical))
        correctedlogical = bits_mask(parities(X, physical))
        @test correctedlogical == xor(rawlogical, recoverylogical)
        @test bits_mask(isodd.(ctx.seed_logicals * UInt8.(seed))) == rawlogical
        push!(histories, (;s0mask, rawlogical, correctedlogical, count=count(seed),
                           physical_signature=bits_mask(parities(H, physical)) |
                               (correctedlogical << n)))
        if recoverylogical != 0 && nonzero_recovery_seed === nothing
            nonzero_recovery_seed = rawmask
        end
        s0mask == 0 && get!(loop_seeds, correctedlogical, rawmask)
    end
    @test Set(keys(loop_seeds)) == Set(0:3)
    @test nonzero_recovery_seed !== nothing

    # Check original raw columns before any same-signature aggregation.
    transfer_signatures = raw_signatures(ctx, ctx.transfer_columns)
    for (local_column, full_column) in enumerate(ctx.transfer_columns)
        @test ctx.transfer_topology.raw_detector_matrix[:,local_column] ==
            ctx.full_topology.raw_detector_matrix[:,full_column]
        @test ctx.transfer_topology.raw_logical_matrix[:,local_column] ==
            ctx.full_topology.raw_logical_matrix[:,full_column]
    end
    for (local_column, full_column) in enumerate(ctx.seed_columns)
        @test ctx.seed_topology.raw_logical_matrix[:,local_column] ==
            ctx.full_topology.raw_logical_matrix[:,full_column]
    end

    # Deterministic physical records include logical loops and cancellation of
    # identical transfer columns; representative posteriors include nonzero s0/sf.
    duplicate_pair = first((i,j) for i in eachindex(transfer_signatures)
        for j in i+1:length(transfer_signatures) if transfer_signatures[i] == transfer_signatures[j])
    fixtures = BitVector[]
    for mask in 0:3
        raw = falses(size(ctx.model.matrix, 2))
        raw[ctx.seed_columns] .= bits_of(loop_seeds[mask], nseed)
        push!(fixtures, raw)
        @test physical_record(ctx, raw).true_mask == mask
    end
    raw = falses(size(ctx.model.matrix, 2))
    raw[ctx.seed_columns] .= bits_of(nonzero_recovery_seed, nseed)
    raw[ctx.transfer_columns[1]] = true
    push!(fixtures, raw)
    canceled = copy(raw)
    canceled[ctx.transfer_columns[collect(duplicate_pair)]] .⊻= true
    @test physical_record(ctx, canceled).final_syndrome == physical_record(ctx, raw).final_syndrome
    @test physical_record(ctx, canceled).true_mask == physical_record(ctx, raw).true_mask

    max_posterior_error = 0.
    for p in (0.02, 0.1)
        decoder = build_decoder(ctx, p)
        @test decoder.seed.settings.max_bond_dimension == 240
        seed_joint = zeros(1 << (n0+2))
        for history in histories
            state = history.s0mask | (history.rawlogical << n0)
            seed_joint[state+1] += p^history.count * (1-p)^(nseed-history.count)
        end
        delta = zeros(1 << (n+2)); delta[1] = 1
        transfer_joint = transfer_dp(delta, transfer_signatures, p)
        @test sum(transfer_joint) ≈ 1 atol=2e-14

        # Direct physical oracle: start with every repaired seed history having
        # the observed s0, then add each of the 104 raw transfer faults by DP.
        for raw in (fixtures[1], fixtures[end])
            record = physical_record(ctx, raw)
            s0, sf = bits_mask(record.seed_syndrome), bits_mask(record.final_syndrome)
            initial = zeros(1 << (n+2))
            for history in histories
                history.s0mask == s0 || continue
                initial[history.physical_signature+1] +=
                    p^history.count * (1-p)^(nseed-history.count)
            end
            direct_joint = transfer_dp(initial, transfer_signatures, p)
            expected_weights = [direct_joint[(sf | (ell << n))+1] for ell in 0:3]
            expected = expected_weights ./ sum(expected_weights)
            observed = (;record.seed_syndrome, record.seed_recovery,
                         record.recovery_mask, record.final_syndrome)
            result = decode_record(decoder, observed)
            @test result.valid
            err = maximum(abs.(result.probabilities - expected))
            max_posterior_error = max(max_posterior_error, err)
            @test result.probabilities ≈ expected atol=2e-9 rtol=2e-8
            @test result.logical_mask == argmax(expected)-1
            @test result.log_joint_syndrome ≈ log(sum(expected_weights)) atol=2e-8
            seed_weights = [seed_joint[(s0 | (ell << n0))+1] for ell in 0:3]
            transfer_weights = [transfer_joint[(sf | (ell << n))+1] for ell in 0:3]
            @test result.seed_result.probabilities ≈ seed_weights ./ sum(seed_weights) atol=2e-10
            @test result.transfer_result.probabilities ≈ transfer_weights ./ sum(transfer_weights) atol=2e-9
            @test result.seed_result.log_syndrome ≈ log(sum(seed_weights)) atol=2e-9
            @test result.transfer_result.log_syndrome ≈ log(sum(transfer_weights)) atol=2e-8
            @test_throws ArgumentError decode_record(decoder,
                merge(observed, (recovery_mask=xor(record.recovery_mask,1),)))
        end

        # Equal RNG seeds reproduce the established sampled-record ordering.
        family = build_toric3d_z_feedback_decoders(ctx.code)
        rng1, rng2 = Xoshiro(2048), Xoshiro(2048)
        for _ in 1:3
            sampled = sample_record(ctx, p, rng1)
            old = sample_toric3d_diagonal_z_faults(rng2, ctx.model.encoder; p)
            flattened = vcat(old.seed, [old.transfer.z[slot, step]
                for (step, op) in enumerate(ctx.model.encoder.operations) for slot in eachindex(op.qubits)])
            @test sampled.raw_faults == flattened
            physical = physical_record(ctx, sampled.raw_faults)
            @test sampled.final_syndrome == physical.final_syndrome
            @test sampled.true_mask == physical.true_mask
            corrected = Toric3DDiagonalZFaultRecord(old.seed .⊻ sampled.seed_recovery, old.transfer)
            replay = replay_toric3d_z_strength(ctx.model.encoder, corrected, family; protocol=:none)
            @test replay.frame.z == physical.frame
            @test replay.syndrome == sampled.final_syndrome
            @test replay.failure == sampled.mwpm_failure
            @test sampled.true_mask == xor(xor(sampled.raw_seed_mask, sampled.recovery_mask), sampled.transfer_mask)
        end
    end
    println("Maximum L=2 posterior absolute error against independent physical DP: ", max_posterior_error)
    @test_throws ArgumentError build_context(1)
    @test_throws ArgumentError build_decoder(ctx, 0.5)
    @test_throws ArgumentError sample_record(ctx, NaN, Xoshiro(1))
end

@testset "repaired-seed experiment API" begin
    for name in (:build_context, :build_decoder, :sample_record, :decode_record,
                 :combine_posteriors)
        @test isdefined(RepairedSeedTN, name)
    end
end
