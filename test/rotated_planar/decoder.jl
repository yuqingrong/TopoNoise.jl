using Test
using PythonCall

@testset "Rotated planar PyMatching backend" begin
    decoder_api = (
        :MatchingDecoders,
        :build_matching_decoders,
        :decode_logical_parities,
        :LogicalFailurePoint,
        :estimate_logical_failure,
    )
    for name in decoder_api
        @test isdefined(TopoNoise, name)
    end

    if all(name -> isdefined(TopoNoise, name), decoder_api)
        sys = pyimport("sys")
        @test !pyconvert(Bool, sys.modules.__contains__("pymatching"))

        for orientation in (:x_ns, :x_ew)
            code = RotatedPlanarCode(3; boundary_orientation=orientation)
            @test maximum(sum(a_s_check_matrix(code); dims=1)) <= 2
            @test maximum(sum(b_p_check_matrix(code); dims=1)) <= 2

            decoders = build_matching_decoders(code)
            @test decoders isa MatchingDecoders
            @test decoders.backend === :PyMatching
            @test v"2.3" <= decoders.backend_version < v"3"
            @test pyconvert(Int, decoders.x_error_matcher.num_detectors) == 4
            @test pyconvert(Int, decoders.x_error_matcher.num_fault_ids) == 1
            @test pyconvert(Int, decoders.z_error_matcher.num_detectors) == 4
            @test pyconvert(Int, decoders.z_error_matcher.num_fault_ids) == 1
        end
        @test pyconvert(Bool, sys.modules.__contains__("pymatching"))

        nongraphlike = RotatedPlanarCode(
            3, :x_ns, [[1], [1], [1]], [[2]])
        @test_throws ArgumentError build_matching_decoders(nongraphlike)
    end
end

function _scalar_failure_counts(rng, encoder, noise, shots)
    code = encoder.code
    decoders = build_matching_decoders(code)
    x_failures = 0
    z_failures = 0
    any_failures = 0
    for _ in 1:shots
        faults = sample_fault_record(rng, encoder, noise)
        frame = propagate_pauli_frame(encoder, faults)
        syndrome = measure_syndrome(code, frame)
        prediction = decode_logical_parities(
            decoders, reshape(syndrome.a_s, 1, :), reshape(syndrome.b_p, 1, :))
        actual_x = isodd(count(@view frame.x[logical_z_support(code)]))
        actual_z = isodd(count(@view frame.z[logical_x_support(code)]))
        failed_x = xor(only(prediction.logical_x), actual_x)
        failed_z = xor(only(prediction.logical_z), actual_z)
        x_failures += failed_x
        z_failures += failed_z
        any_failures += failed_x | failed_z
    end
    return (x=x_failures, z=z_failures, any=any_failures)
end

@testset "Chunked logical-failure estimation" begin
    @testset "validation precedes simulation" begin
        encoder = rotated_planar_encoder(RotatedPlanarCode(3))
        noise = CircuitPauliNoise(0)
        @test_throws ArgumentError estimate_logical_failure(
            MersenneTwister(1), encoder, noise; shots=0)
        @test_throws ArgumentError estimate_logical_failure(
            MersenneTwister(1), encoder, noise; shots=-1)
        @test_throws ArgumentError estimate_logical_failure(
            MersenneTwister(1), encoder, noise; shots=1, batch_size=0)
        @test_throws ArgumentError estimate_logical_failure(
            MersenneTwister(1), encoder, noise; shots=1, batch_size=-1)
        @test_throws ArgumentError estimate_logical_failure(
            MersenneTwister(1), encoder, noise; shots=1, seed=1.5)
        @test_throws ArgumentError estimate_logical_failure(
            MersenneTwister(1), encoder, noise; shots=1, seed="not-an-integer")
    end

    @testset "p=0 has exact zero failures across all metadata axes" begin
        for orientation in (:x_ns, :x_ew), construction in (:as, :bp),
            logical_state in (:zero, :one, :plus, :minus),
            clock in (:gate_layer, :plaquette)
            code = RotatedPlanarCode(3; boundary_orientation=orientation)
            encoder = rotated_planar_encoder(
                code; construction=construction, logical_state=logical_state)
            point = estimate_logical_failure(
                MersenneTwister(80), encoder, CircuitPauliNoise(0; clock=clock);
                shots=3, batch_size=2, seed=UInt32(80))

            @test point isa LogicalFailurePoint
            @test point.distance == 3
            @test point.boundary_orientation === orientation
            @test point.construction === construction
            @test point.logical_state === logical_state
            @test point.clock === clock
            @test point.p_x == 0.0
            @test point.p_z == 0.0
            @test point.shots == 3
            @test point.seed == 80
            @test point.logical_x_failures == 0
            @test point.logical_z_failures == 0
            @test point.any_logical_failures == 0
            @test point.state_failures == 0
            @test point.logical_x_failure_rate == 0.0
            @test point.logical_z_failure_rate == 0.0
            @test point.any_logical_failure_rate == 0.0
            @test point.state_failure_rate == 0.0
            @test point.logical_x_standard_error == 0.0
            @test point.logical_z_standard_error == 0.0
            @test point.any_logical_standard_error == 0.0
            @test point.state_standard_error == 0.0
        end
    end

    @testset "vector engine agrees with the scalar frame oracle" begin
        code = RotatedPlanarCode(3; boundary_orientation=:x_ew)
        encoder = rotated_planar_encoder(
            code; construction=:as, logical_state=:plus)
        noise = CircuitPauliNoise(0; p_x=0.08, p_z=0.31, clock=:plaquette)
        shots = 17
        expected = _scalar_failure_counts(
            MersenneTwister(902), encoder, noise, shots)
        actual = estimate_logical_failure(
            MersenneTwister(902), encoder, noise;
            shots=shots, batch_size=1, seed=nothing)

        @test actual.logical_x_failures == expected.x
        @test actual.logical_z_failures == expected.z
        @test actual.any_logical_failures == expected.any
        @test actual.state_failures == expected.z
        @test actual.seed === nothing
    end

    @testset "partial batches, RNG ownership, rates, and state selection" begin
        code = RotatedPlanarCode(3)
        noise = CircuitPauliNoise(0; p_x=0.04, p_z=0.37, clock=:gate_layer)
        plus_encoder = rotated_planar_encoder(
            code; construction=:bp, logical_state=:plus)
        first = estimate_logical_failure(
            MersenneTwister(771), plus_encoder, noise;
            shots=23, batch_size=7, seed=771)
        repeated = estimate_logical_failure(
            MersenneTwister(771), plus_encoder, noise;
            shots=23, batch_size=7, seed=771)
        relabeled = estimate_logical_failure(
            MersenneTwister(771), plus_encoder, noise;
            shots=23, batch_size=7, seed=999)

        for field in fieldnames(LogicalFailurePoint)
            @test getfield(first, field) == getfield(repeated, field)
            field === :seed || @test getfield(first, field) == getfield(relabeled, field)
        end
        @test relabeled.seed == 999
        @test first.p_x == 0.04
        @test first.p_z == 0.37
        @test 0 <= first.logical_x_failures <= first.shots
        @test 0 <= first.logical_z_failures <= first.shots
        @test max(first.logical_x_failures, first.logical_z_failures) <=
              first.any_logical_failures <= first.shots
        @test first.logical_x_failure_rate == first.logical_x_failures / first.shots
        @test first.logical_z_failure_rate == first.logical_z_failures / first.shots
        @test first.any_logical_failure_rate == first.any_logical_failures / first.shots
        @test first.state_failure_rate == first.state_failures / first.shots
        @test first.logical_x_standard_error ≈ sqrt(
            first.logical_x_failure_rate * (1 - first.logical_x_failure_rate) /
            first.shots)
        @test first.logical_z_standard_error ≈ sqrt(
            first.logical_z_failure_rate * (1 - first.logical_z_failure_rate) /
            first.shots)
        @test first.any_logical_standard_error ≈ sqrt(
            first.any_logical_failure_rate * (1 - first.any_logical_failure_rate) /
            first.shots)
        @test first.state_standard_error ≈ sqrt(
            first.state_failure_rate * (1 - first.state_failure_rate) /
            first.shots)
        @test first.logical_x_failures != first.logical_z_failures
        @test first.state_failures == first.logical_z_failures

        zero_encoder = rotated_planar_encoder(
            code; construction=:bp, logical_state=:zero)
        zero_point = estimate_logical_failure(
            MersenneTwister(772), zero_encoder, noise;
            shots=23, batch_size=7)
        @test zero_point.logical_x_failures != zero_point.logical_z_failures
        @test zero_point.state_failures == zero_point.logical_x_failures

        @test TopoNoise._binomial_standard_error(0.0, 7) == 0.0
        @test TopoNoise._binomial_standard_error(1.0, 7) == 0.0
    end

    @testset "distance five exercises the production binary path" begin
        code = RotatedPlanarCode(5; boundary_orientation=:x_ew)
        encoder = rotated_planar_encoder(
            code; construction=:as, logical_state=:minus)
        point = estimate_logical_failure(
            MersenneTwister(5), encoder,
            CircuitPauliNoise(0; clock=:plaquette);
            shots=1, batch_size=1)
        @test point.shots == 1
        @test point.any_logical_failures == 0
    end
end

@testset "Syndrome-only batch logical decoding" begin
    a_s_fixtures = Dict(
        :x_ns => Bool[
            1 1 0 1 1 0 0 0 0
            0 0 0 0 1 1 0 1 1
            0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0
        ],
        :x_ew => Bool[
            0 1 1 0 1 1 0 0 0
            0 0 0 1 1 0 1 1 0
            0 0 0 0 0 1 0 0 1
            1 0 0 1 0 0 0 0 0
        ],
    )
    b_p_fixtures = Dict(
        :x_ns => Bool[
            0 1 1 0 1 1 0 0 0
            0 0 0 1 1 0 1 1 0
            1 0 0 1 0 0 0 0 0
            0 0 0 0 0 1 0 0 1
        ],
        :x_ew => Bool[
            0 0 0 0 1 1 0 1 1
            1 1 0 1 1 0 0 0 0
            0 1 1 0 0 0 0 0 0
            0 0 0 0 0 0 1 1 0
        ],
    )
    logical_z_fixtures = Dict(:x_ns => [1, 4, 7], :x_ew => [1, 2, 3])
    logical_x_fixtures = Dict(:x_ns => [1, 2, 3], :x_ew => [3, 6, 9])

    for orientation in (:x_ns, :x_ew)
        code = RotatedPlanarCode(3; boundary_orientation=orientation)
        decoders = build_matching_decoders(code)
        a_s = permutedims(a_s_fixtures[orientation])
        b_p = permutedims(b_p_fixtures[orientation])
        zero_a_s = falses(9, 4)
        zero_b_p = falses(9, 4)
        expected_x = BitVector(qubit in logical_z_fixtures[orientation] for qubit in 1:9)
        expected_z = BitVector(qubit in logical_x_fixtures[orientation] for qubit in 1:9)

        x_only = decode_logical_parities(decoders, a_s, zero_b_p)
        @test x_only.logical_x == expected_x
        @test x_only.logical_z == falses(9)

        z_only = decode_logical_parities(decoders, zero_a_s, b_p)
        @test z_only.logical_x == falses(9)
        @test z_only.logical_z == expected_z

        y = decode_logical_parities(decoders, a_s, b_p)
        @test y.logical_x == expected_x
        @test y.logical_z == expected_z

        # Minimum logical strings are undetectable. Their actual parities are
        # one, while a final-syndrome-only decoder necessarily predicts zero.
        logical_x_frame = falses(9)
        logical_x_frame[logical_x_fixtures[orientation]] .= true
        logical_z_frame = falses(9)
        logical_z_frame[logical_z_fixtures[orientation]] .= true
        logical_x_syndrome = measure_syndrome(
            code, PauliFrame(logical_x_frame, falses(9)))
        logical_z_syndrome = measure_syndrome(
            code, PauliFrame(falses(9), logical_z_frame))
        @test !any(logical_x_syndrome.a_s)
        @test !any(logical_z_syndrome.b_p)
        @test isodd(count(logical_x_frame[logical_z_fixtures[orientation]]))
        @test isodd(count(logical_z_frame[logical_x_fixtures[orientation]]))
        undetectable = decode_logical_parities(
            decoders,
            reshape(logical_x_syndrome.a_s, 1, :),
            reshape(logical_z_syndrome.b_p, 1, :),
        )
        @test !only(undetectable.logical_x)
        @test !only(undetectable.logical_z)
        @test xor(only(undetectable.logical_x), true)
        @test xor(only(undetectable.logical_z), true)

        @test_throws DimensionMismatch decode_logical_parities(
            decoders, falses(2, 3), falses(2, 4))
        @test_throws DimensionMismatch decode_logical_parities(
            decoders, falses(2, 4), falses(1, 4))
    end

    # The public decoding boundary consists solely of the decoder handle and
    # the two final CSS syndrome matrices; actual errors/parities stay in Julia.
    decoder_methods = collect(methods(decode_logical_parities))
    @test length(decoder_methods) == 1
    @test only(decoder_methods).nargs == 4
    @test hasmethod(
        decode_logical_parities,
        Tuple{MatchingDecoders, BitMatrix, BitMatrix})
end
