using TopoNoise, Random, Statistics, Printf, SHA

row_reference(matrix, frame) = BitVector(isodd(count((@view matrix[row, :]) .& frame))
    for row in axes(matrix, 1))
function hot_time(f, matrix, frame)
    f(matrix, frame)
    median([@elapsed f(matrix, frame) for _ in 1:7])
end
function paired(encoder, faults, family; record_history=true)
    original = deepcopy(faults)
    results = (final_only=replay_toric3d_h6_full_prefix_x(encoder, faults, family;
        feedback=false, record_history),
        full_prefix=replay_toric3d_h6_full_prefix_x(encoder, faults, family;
        feedback=true, record_history))
    @assert faults == original
    results
end
function main()
    path = joinpath(dirname(pathof(TopoNoise)), "toric3d", "decoder.jl")
    initial_hash = bytes2hex(sha256(read(path)))
    println("OPTIMIZED_SOURCE_SHA256=", initial_hash)
    println("MICROBENCH: actual L11 plaquette matrices; median 7 hot calls")
    println("prefix,rows,columns,density,weight,equal,row_ms,packed_ms,speedup,row_bytes,packed_bytes")
    rng = MersenneTwister(20261015)
    equality_cases = 0
    for r in 1:11
        matrix = BitMatrix(plaquette_check_matrix(Toric3DEncodedSlice(r, 11, 11)))
        for density in (0.0, 0.01, 0.05, 0.1, 0.5, 1.0), repetition in 1:3
            frame = BitVector(rand(rng, size(matrix, 2)) .< density)
            @assert row_reference(matrix, frame) == TopoNoise._toric3d_matrix_syndrome(matrix, frame)
            equality_cases += 1
            if r in (1, 5, 11) && repetition == 1 && density != .1
                old = hot_time(row_reference, matrix, frame)
                new = hot_time(TopoNoise._toric3d_matrix_syndrome, matrix, frame)
                oldbytes = @allocated row_reference(matrix, frame)
                newbytes = @allocated TopoNoise._toric3d_matrix_syndrome(matrix, frame)
                @printf("%d,%d,%d,%.2f,%d,true,%.6f,%.6f,%.2f,%d,%d\n",
                    r, size(matrix, 1), size(matrix, 2), density, count(frame),
                    1000old, 1000new, old/new, oldbytes, newbytes)
                flush(stdout)
            end
        end
    end
    println("ALL_PREFIX_RANDOM_EQUALITY_CASES=", equality_cases)
    code = Toric3DEncodedSlice(11, 11, 11)
    encoder = toric3d_encoded_encoder(code)
    family = build_toric3d_x_feedback_decoders(code; p=.05)
    faults = [sample_toric3d_faults(MersenneTwister(seed), encoder; p_x=.05)
              for seed in (314159, 271828)]
    paired(encoder, faults[1], family)
    optimized = Any[]; optimized_seconds = Float64[]
    for record in faults
        seconds = @elapsed result = paired(encoder, record, family)
        push!(optimized, result); push!(optimized_seconds, seconds)
    end
    optimized_nohistory = @elapsed paired(encoder, faults[1], family; record_history=false)
    println("OPTIMIZED_HOT_PAIR_SECONDS=", optimized_seconds)
    @printf("OPTIMIZED_NOHISTORY_PAIR_SECONDS=%.6f\n", optimized_nohistory)
    flush(stdout)
    # Restore the previous row function only inside this verification process.
    # Production files and other processes are never changed by this swap.
    @eval TopoNoise function _toric3d_matrix_syndrome(matrix::BitMatrix, frame::BitVector)
        BitVector(isodd(count((@view matrix[row, :]) .& frame)) for row in axes(matrix, 1))
    end
    Base.invokelatest(paired, encoder, faults[1], family)
    previous = Any[]; previous_seconds = Float64[]
    for record in faults
        seconds = @elapsed result = Base.invokelatest(paired, encoder, record, family)
        push!(previous, result); push!(previous_seconds, seconds)
    end
    previous_nohistory = @elapsed Base.invokelatest(paired, encoder, faults[1], family; record_history=false)
    @assert optimized == previous
    @assert bytes2hex(sha256(read(path))) == initial_hash
    println("PREVIOUS_HOT_PAIR_SECONDS=", previous_seconds)
    @printf("PREVIOUS_NOHISTORY_PAIR_SECONDS=%.6f\n", previous_nohistory)
    println("EXACT_ALL_RESULTS_AND_HISTORIES_EQUAL=true")
    println("SOURCE_HASH_UNCHANGED_DURING_VERIFICATION=true")
    for i in eachindex(faults)
        @printf("PAIR_%d_SPEEDUP=%.3f\n", i, previous_seconds[i]/optimized_seconds[i])
        println("PAIR_", i, "_FAILURES=", (previous[i].final_only.failure, previous[i].full_prefix.failure))
    end
    @printf("NOHISTORY_SPEEDUP=%.3f\n", previous_nohistory/optimized_nohistory)
    println("COMPLETE")
end
main()
