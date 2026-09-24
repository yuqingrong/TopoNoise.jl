using Test, Random
import Yao
using LinearAlgebra: dot

@testset "Post-encoding IID Pauli noise" begin
    # A rejected clock must fail this test before the new sampling branch exists.
    candidate = try
        CircuitPauliNoise(0.13; p_x=0.13, p_z=0, clock=:post_encoding)
    catch error
        error
    end
    @test candidate isa CircuitPauliNoise

    if candidate isa CircuitPauliNoise
        noise = candidate
        @testset "one Bernoulli trial per data qubit, independent of encoder" begin
            for d in (3, 5), orientation in (:x_ns, :x_ew)
                code = RotatedPlanarCode(d; boundary_orientation=orientation)
                expected_x = BitMatrix(rand(MersenneTwister(812), 41, d^2) .< 0.13)
                records = CircuitFaultRecord[]
                for construction in (:as, :bp)
                    encoder = rotated_planar_encoder(code; construction)
                    @test verify_encoder_tableau(encoder)
                    x, z = TopoNoise._sample_frame_batch(
                        MersenneTwister(812), encoder, noise, 41)
                    @test x == expected_x
                    @test z == falses(41, d^2)
                    record = sample_fault_record(MersenneTwister(91), encoder, noise)
                    push!(records, record)
                    @test length(record.steps) == 1
                    @test only(record.steps).support == Tuple(1:d^2)
                    frame = propagate_pauli_frame(encoder, record)
                    @test frame.x == collect(only(record.steps).x)
                    @test frame.z == falses(d^2)
                end
                @test first(records).steps == last(records).steps
            end
        end

        @testset "final single faults never propagate back through the encoder" begin
            for orientation in (:x_ns, :x_ew), construction in (:as, :bp)
                code = RotatedPlanarCode(3; boundary_orientation=orientation)
                encoder = rotated_planar_encoder(code; construction)
                clean = sample_fault_record(MersenneTwister(8), encoder,
                    CircuitPauliNoise(0; clock=:post_encoding))
                for q in 1:9, pauli in (:X, :Z, :Y)
                    record = with_pauli_fault(clean, 1, q, pauli)
                    expected_x = falses(9)
                    expected_z = falses(9)
                    expected_x[q] = pauli in (:X, :Y)
                    expected_z[q] = pauli in (:Z, :Y)
                    expected = PauliFrame(expected_x, expected_z)
                    @test propagate_pauli_frame(encoder, record) == expected
                    @test sample_yao_syndrome(MersenneTwister(9), code, encoder, record) ==
                          measure_syndrome(code, expected)
                end
                @test_throws ArgumentError with_pauli_fault(clean, 2, 1, :X)
                @test_throws ArgumentError with_pauli_fault(clean, 1, 10, :X)
            end
        end

        @testset "As and Bp give the same noisy Yao state after ideal preparation" begin
            for orientation in (:x_ns, :x_ew), state in (:zero, :one, :plus, :minus)
                code = RotatedPlanarCode(3; boundary_orientation=orientation)
                registers = []
                for construction in (:as, :bp)
                    encoder = rotated_planar_encoder(code; construction, logical_state=state)
                    record = sample_fault_record(MersenneTwister(404), encoder,
                        CircuitPauliNoise(0.17; p_z=0.23, clock=:post_encoding))
                    actual = Yao.zero_state(9)
                    TopoNoise._apply_yao_encoder_and_faults!(actual, 9, encoder, record)
                    expected = Yao.zero_state(9)
                    Yao.apply!(expected, yao_encoder(encoder))
                    for q in 1:9
                        only(record.steps).x[q] && Yao.apply!(expected, Yao.put(9, q => Yao.X))
                        only(record.steps).z[q] && Yao.apply!(expected, Yao.put(9, q => Yao.Z))
                    end
                    @test abs(dot(Yao.statevec(expected), Yao.statevec(actual))) ≈ 1 atol=1e-10
                    push!(registers, actual)
                end
                @test abs(dot(Yao.statevec(first(registers)),
                              Yao.statevec(last(registers)))) ≈ 1 atol=1e-10
            end
        end

        @testset "batched decoding agrees across constructions and preserves counts" begin
            code = RotatedPlanarCode(3)
            points = LogicalFailurePoint[]
            for construction in (:as, :bp)
                encoder = rotated_planar_encoder(code; construction, logical_state=:zero)
                clean = estimate_logical_failure(MersenneTwister(21), encoder,
                    CircuitPauliNoise(0; clock=:post_encoding); shots=17, batch_size=7)
                @test clean.any_logical_failures == 0
                point = estimate_logical_failure(MersenneTwister(22), encoder,
                    noise; shots=101, batch_size=23, seed=22)
                push!(points, point)
                @test point.clock === :post_encoding
                @test point.shots == 101
                @test point.logical_z_failures == 0
                @test point.state_failures == point.logical_x_failures == point.any_logical_failures
                @test point.logical_x_failure_rate == point.logical_x_failures / 101
                @test point.logical_x_standard_error ≈ sqrt(
                    point.logical_x_failure_rate * (1 - point.logical_x_failure_rate) / 101)
            end
            @test first(points).logical_x_failures == last(points).logical_x_failures
            @test 0 < first(points).logical_x_failures < 101

            encoder = rotated_planar_encoder(code)
            decoders = build_matching_decoders(code)
            for q in 1:9
                x = falses(9)
                x[q] = true
                syndrome = measure_syndrome(code, PauliFrame(x, falses(9)))
                predicted = decode_logical_parities(decoders,
                    reshape(syndrome.a_s, 1, :), reshape(syndrome.b_p, 1, :))
                @test only(predicted.logical_x) == (q in logical_z_support(code))
                @test !only(predicted.logical_z)
            end
            for p in (-0.1, 0.5, Inf, NaN)
                @test_throws ArgumentError CircuitPauliNoise(p; clock=:post_encoding)
            end
            broken = CircuitFaultRecord(:post_encoding, 3, :x_ns, :bp, :zero,
                [CircuitFaultStep(1, :post_encoding, 0, [1], [false], [false])])
            @test_throws ArgumentError propagate_pauli_frame(encoder, broken)
        end
    end
end

if !isdefined(@__MODULE__, :RotatedPlanarConstructionComparison)
    include(joinpath(@__DIR__, "..", "..", "examples", "compare_rotated_planar_constructions.jl"))
end

@testset "Comparison CLI accepts the post-encoding model" begin
    settings = try
        RotatedPlanarConstructionComparison._parse_arguments([
            "--clock", "post_encoding", "--distances", "3", "--shots", "2"])
    catch error
        error
    end
    @test settings isa Dict
    if settings isa Dict
        @test settings[:clock] === :post_encoding
        @test settings[:shots] == 2
        @test settings[:distances] == [3]
    end
end
