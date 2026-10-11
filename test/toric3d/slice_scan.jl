using Test
using Random
import TOML

@testset "Slice-transfer comparison API" begin
    @test isdefined(TopoNoise, :compare_toric3d_slice_feedback)
    @test isdefined(TopoNoise, :save_toric3d_slice_comparison)
    @test isdefined(TopoNoise, :toric3d_slice_diagnostics)
end

if isdefined(TopoNoise, :compare_toric3d_slice_feedback)
    include(joinpath(@__DIR__,"..","..","examples","compare_toric3d_slice_feedback.jl"))
    @testset "Slice comparison example arguments" begin
        output,errors = IOBuffer(),IOBuffer()
        @test Toric3DSliceFeedbackExample.main(["--help"];io=output,error_io=errors) == 0
        @test occursin("No final global recovery",String(take!(output)))
        @test isempty(String(take!(errors)))
        for arguments in (["--shots","0"],["--routing","other"],["--sizes","2,2"],
                ["--error-rates",".03,.01"],["--seed","-1"],["--unknown","1"])
            @test Toric3DSliceFeedbackExample.main(arguments;io=output,error_io=errors) == 1
            @test !isempty(String(take!(errors)))
        end
        options = Toric3DSliceFeedbackExample._parse_arguments(
            ["--sizes","2,3","--routing","full_path","--shots","12"])
        @test options[:sizes] == [2,3]
        @test options[:routing] === :full_path
        @test options[:shots] == 12
    end
end

if isdefined(TopoNoise, :compare_toric3d_slice_feedback)
    @testset "Slice statistics use independent shots and paired differences" begin
        m = TopoNoise._toric3d_slice_metric([0., .5, 1.], [.5, .5, .5])
        @test m.unassisted_mean == m.slice_feedback_mean == .5
        @test m.unassisted_standard_error ≈ sqrt(1 / 12)
        @test m.slice_feedback_standard_error == 0
        @test m.difference == 0
        @test m.difference_standard_error ≈ sqrt(1 / 12)
        @test m.difference_ci[1] < 0 < m.difference_ci[2]
    end

    @testset "Plaquette clusters connect through shared cubes" begin
        code = Toric3DSliceLattice(3, 2, 2)
        # Face order is t,u,v,pair. B12(t,0,0) touches the same cube as
        # B12(t-1,0,0), because s3=(0,0).
        faces = falses(36)
        @test TopoNoise._toric3d_slice_cluster_extent(code, faces) == 0
        faces[[1, 13, 25]] .= true
        @test TopoNoise._toric3d_slice_cluster_extent(code, faces) == 1
        faces[13] = false
        @test TopoNoise._toric3d_slice_cluster_extent(code, faces) == 1 / 3
        faces[13] = true
        faces[25] = false
        @test TopoNoise._toric3d_slice_cluster_extent(code, faces) == 2 / 3
        @test_throws DimensionMismatch TopoNoise._toric3d_slice_cluster_extent(code, falses(35))

        rectangular = Toric3DSliceLattice(2,3,3)
        shifted = falses(54)
        shifted[[8,29]] .= true # B13 below at (0,2), above at (0,0).
        @test TopoNoise._toric3d_slice_cluster_extent(rectangular,shifted) == 1
        shifted .= false
        shifted[[21,30]] .= true # B23 below at (2,0), above at (0,0).
        @test TopoNoise._toric3d_slice_cluster_extent(rectangular,shifted) == 1
        shifted .= false
        shifted[[2,29]] .= true # Unshifted B13 faces share no cube.
        @test TopoNoise._toric3d_slice_cluster_extent(rectangular,shifted) == .5

        syndromes = (stars=BitVector([1,0,0,0,0,0,0,0,1,1,1,1]),
            plaquettes=falses(36), boundary=falses(12), frontier=BitVector([1,0,0,0]))
        metrics = toric3d_slice_diagnostics(code, syndromes)
        @test metrics.bulk_star_density == 1 / 8
        @test metrics.plaquette_density == 0
        @test metrics.frontier_density == 1 / 4
        @test metrics.plaquette_cluster_extent == 0
        thin = Toric3DSliceLattice(1, 2, 2)
        @test toric3d_slice_diagnostics(thin,
            (stars=trues(4), plaquettes=falses(12), boundary=falses(12), frontier=falses(4))).bulk_star_density == 0
    end

    @testset "Slice scan preserves zero noise, pairing, and reproducibility" begin
        args = (sizes=[(2,2,2)], error_rates=[0., .03], shots=8, seed=71)
        a = compare_toric3d_slice_feedback(MersenneTwister(71); args...)
        b = compare_toric3d_slice_feedback(MersenneTwister(71); args...)
        @test length(a.points) == 4
        @test [p.metrics for p in a.points] == [p.metrics for p in b.points]
        for point in a.points
            @test point.shots == 8
            for metric in values(point.metrics)
                @test 0 <= metric.unassisted_mean <= 1
                @test 0 <= metric.slice_feedback_mean <= 1
                @test metric.difference ≈ metric.slice_feedback_mean - metric.unassisted_mean atol=1e-15
            end
            if point.p == 0
                @test all(m.unassisted_mean == m.slice_feedback_mean == 0 for m in values(point.metrics))
            end
            if point.channel === :z_only
                @test point.metrics.bulk_star_density.difference == 0
                @test point.metrics.bulk_star_density.difference_standard_error == 0
                @test point.metrics.plaquette_density.unassisted_mean == point.metrics.plaquette_density.slice_feedback_mean == 0
                @test point.metrics.frontier_density.unassisted_mean == point.metrics.frontier_density.slice_feedback_mean == 0
            end
        end
        @test_throws ArgumentError compare_toric3d_slice_feedback(MersenneTwister(1); sizes=[], shots=1)
        @test_throws ArgumentError compare_toric3d_slice_feedback(MersenneTwister(1); shots=0)
        @test_throws ArgumentError compare_toric3d_slice_feedback(MersenneTwister(1); error_rates=[.03,.01])
        @test_throws ArgumentError compare_toric3d_slice_feedback(MersenneTwister(1); routing=:unknown)
        full = compare_toric3d_slice_feedback(MersenneTwister(4);
            sizes=[(1,2,3)],error_rates=[.03],shots=2,routing=:full_path)
        @test full.routing === :full_path
        @test full.dimensions == [(1,2,3)]
        @test all(p.metrics.bulk_star_density.unassisted_mean == p.metrics.bulk_star_density.slice_feedback_mean == 0 for p in full.points)
        mktempdir() do directory
            c = compare_toric3d_slice_feedback(MersenneTwister(71); args..., output_dir=directory)
            raw = joinpath(directory,"toric3d-slice-feedback-raw.csv")
            @test length(readlines(raw)) == 1 + 4length(c.points)
            metadata = TOML.parsefile(joinpath(directory,"toric3d-slice-feedback-metadata.toml"))
            @test metadata["complete"]
            @test metadata["completed_channel_points"] == 4
            @test metadata["seed_metadata"] == 71
            @test metadata["routing"] == "one_step"
            paths = save_toric3d_slice_comparison(c, directory)
            for key in (:csv,:metadata,:svg,:pdf,:png)
                @test isfile(getproperty(paths,key))
                @test filesize(getproperty(paths,key)) > 0
            end
        end
    end
end
