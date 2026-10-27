using TopoNoise, Random, Statistics, Printf, SHA

function column_syndrome(matrix::BitMatrix, frame::BitVector)
    result = falses(size(matrix, 1))
    for j in findall(frame)
        result .⊻= view(matrix, :, j)
    end
    return result
end

function hot_time(f, matrix, frame; samples=7)
    f(matrix, frame)
    times = Float64[]
    for _ in 1:samples
        push!(times, @elapsed f(matrix, frame))
    end
    return median(times)
end

function paired(encoder, faults, family; record_history=true)
    original = deepcopy(faults)
    result = (final_only=replay_toric3d_h6_full_prefix_x(encoder, faults, family;
        feedback=false, record_history),
        full_prefix=replay_toric3d_h6_full_prefix_x(encoder, faults, family;
        feedback=true, record_history))
    @assert faults == original
    return result
end

function main()
    path = joinpath(dirname(pathof(TopoNoise)), "toric3d", "decoder.jl")
    initial_hash = bytes2hex(sha256(read(path)))
    println("SOURCE_SHA256=", initial_hash)
    println("MICROBENCH: median 7 hot calls; single Julia thread; actual L11 plaquette matrices")
    println("prefix,rows,columns,density,weight,equal,row_ms,column_ms,speedup,row_bytes,column_bytes")
    rng = MersenneTwister(20261015)
    micro = NamedTuple[]
    for r in (1, 5, 11)
        code = Toric3DEncodedSlice(r, 11, 11)
        matrix = BitMatrix(plaquette_check_matrix(code))
        for density in (0.0, 0.01, 0.05, 0.5, 1.0)
            frame = BitVector(rand(rng, size(matrix, 2)) .< density)
            expected = TopoNoise._toric3d_matrix_syndrome(matrix, frame)
            actual = column_syndrome(matrix, frame)
            @assert actual == expected
            old = hot_time(TopoNoise._toric3d_matrix_syndrome, matrix, frame)
            new = hot_time(column_syndrome, matrix, frame)
            oldbytes = @allocated TopoNoise._toric3d_matrix_syndrome(matrix, frame)
            newbytes = @allocated column_syndrome(matrix, frame)
            @printf("%d,%d,%d,%.2f,%d,true,%.6f,%.6f,%.2f,%d,%d\n",
                r, size(matrix, 1), size(matrix, 2), density, count(frame),
                1000old, 1000new, old/new, oldbytes, newbytes)
            push!(micro, (r=r, density=density, old=old, new=new))
            flush(stdout)
        end
    end
    equal_cases = 0
    for r in 1:11
        matrix = BitMatrix(plaquette_check_matrix(Toric3DEncodedSlice(r, 11, 11)))
        for density in (0.0, 0.01, 0.05, 0.1, 0.5, 1.0), _ in 1:3
            frame = BitVector(rand(rng, size(matrix, 2)) .< density)
            @assert TopoNoise._toric3d_matrix_syndrome(matrix, frame) == column_syndrome(matrix, frame)
            equal_cases += 1
        end
    end
    println("ALL_PREFIX_RANDOM_EQUALITY_CASES=", equal_cases)
    flush(stdout)

    println("Building existing L11 decoder family, p=0.05, unchanged default settings")
    code = Toric3DEncodedSlice(11, 11, 11)
    encoder = toric3d_encoded_encoder(code)
    build_seconds = @elapsed family = build_toric3d_x_feedback_decoders(code; p=.05)
    @printf("DECODER_BUILD_SECONDS=%.3f\n", build_seconds)
    faults = [sample_toric3d_faults(MersenneTwister(seed), encoder; p_x=.05)
              for seed in (314159, 271828)]
    paired(encoder, faults[1], family)
    old_results = Any[]
    old_times = Float64[]
    for record in faults
        elapsed = @elapsed result = paired(encoder, record, family)
        push!(old_results, result); push!(old_times, elapsed)
    end
    old_nohistory = @elapsed paired(encoder, faults[1], family; record_history=false)
    println("OLD_HOT_PAIR_SECONDS=", old_times)
    @printf("OLD_NOHISTORY_PAIR_SECONDS=%.6f\n", old_nohistory)
    flush(stdout)

    # This replacement exists only in this isolated Julia process. It does not
    # write package files or change any other running production process.
    @eval TopoNoise function _toric3d_matrix_syndrome(matrix::BitMatrix, frame::BitVector)
        result = falses(size(matrix, 1))
        for j in findall(frame)
            result .⊻= view(matrix, :, j)
        end
        return result
    end
    Base.invokelatest(paired, encoder, faults[1], family)
    new_results = Any[]
    new_times = Float64[]
    for record in faults
        elapsed = @elapsed result = Base.invokelatest(paired, encoder, record, family)
        push!(new_results, result); push!(new_times, elapsed)
    end
    new_nohistory = @elapsed Base.invokelatest(paired, encoder, faults[1], family; record_history=false)
    @assert old_results == new_results
    @assert bytes2hex(sha256(read(path))) == initial_hash
    println("NEW_HOT_PAIR_SECONDS=", new_times)
    @printf("NEW_NOHISTORY_PAIR_SECONDS=%.6f\n", new_nohistory)
    println("EXACT_ALL_RESULTS_AND_HISTORIES_EQUAL=true")
    println("SOURCE_HASH_UNCHANGED=true")
    for i in eachindex(faults)
        @printf("PAIR_%d_SPEEDUP=%.3f\n", i, old_times[i]/new_times[i])
        println("PAIR_", i, "_FAILURES=", (old_results[i].final_only.failure, old_results[i].full_prefix.failure))
    end
    @printf("NOHISTORY_SPEEDUP=%.3f\n", old_nohistory/new_nohistory)
    println("COMPLETE")
end

main()
