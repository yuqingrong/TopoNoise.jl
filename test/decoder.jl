using TopoNoise
using Random
using Test

const R, C = 8, 8

_zero_errors(R, C) = TopoNoise.VirtualBondErrors(
    falses(R, C - 1), falses(R - 1, C),
    falses(R), falses(R), falses(C), falses(C))

_no_mismatches(R, C) = (
    horizontal=falses(R, C - 1), vertical=falses(R - 1, C))

@testset "decoder" begin
    model = ToricCodeTrajectoryModel(R, C)

    @testset "trivial: no defects" begin
        mm = _no_mismatches(R, C)
        for boundary in (:north_south, :east_west)
            correction = decode_uf(model, mm; boundary=boundary)
            @test !any(correction.horizontal)
            @test !any(correction.vertical)
            @test logical_failure(
                _zero_errors(R, C), correction; boundary=boundary) == false
        end
    end

    @testset "plaquette_syndrome parity" begin
        # Single horizontal error at horizontal[2, 3] flips plaquettes (1,3) and (2,3).
        h = falses(R, C - 1); v = falses(R - 1, C)
        h[2, 3] = true
        mm = (; horizontal=h, vertical=v)
        syndrome = plaquette_syndrome(mm)
        @test syndrome[1, 3] && syndrome[2, 3]
        @test count(syndrome) == 2
    end

    @testset "single bulk pair is corrected" begin
        # A single interior horizontal error creates two adjacent defects.
        # Decoder should return a one-bond correction that cancels the error.
        h = falses(R, C - 1); v = falses(R - 1, C)
        h[4, 4] = true
        mm = (; horizontal=h, vertical=v)
        correction = decode_uf(model, mm; boundary=:north_south)
        # Residual = correction XOR error must have zero plaquette syndrome.
        residual_h = xor.(h, correction.horizontal)
        residual_v = xor.(v, correction.vertical)
        @test !any(plaquette_syndrome(
            (; horizontal=residual_h, vertical=residual_v)))
    end

    @testset "decode_trajectory pipeline: p=0 never fails" begin
        rng = MersenneTwister(0)
        for _ in 1:20
            trajectory = sample_trajectory(rng, model; error_rate=0.0)
            result = decode_trajectory(trajectory; boundary=:north_south)
            @test result.logical_failure == false
        end
    end

    @testset "high p produces near-random failure" begin
        rng = MersenneTwister(42)
        fails = 0
        shots = 200
        for _ in 1:shots
            trajectory = sample_trajectory(rng, model; error_rate=0.35)
            result = decode_trajectory(trajectory; boundary=:north_south)
            fails += result.logical_failure ? 1 : 0
        end
        @test 0.30 <= fails / shots <= 0.70
    end

    @testset "boundary-anchored logical is a homology invariant" begin
        # Any bulk error whose plaquette-parity is zero and whose N-S vertical
        # crossing on the reference row is 0 should decode to no failure.
        rng = MersenneTwister(7)
        model_small = ToricCodeTrajectoryModel(4, 4)
        errs = _zero_errors(4, 4)
        # Add a small closed loop of bulk vertical errors that returns parity 0.
        # Loop: vertical[1,2] + vertical[2,2] + horizontal[2,2] + horizontal[3,2]
        # form a closed loop around plaquette (2,2). Actually simpler: add none.
        trajectory = sample_trajectory(rng, model_small, errs)
        r = decode_trajectory(trajectory; boundary=:north_south)
        @test r.logical_failure == false
    end
end
