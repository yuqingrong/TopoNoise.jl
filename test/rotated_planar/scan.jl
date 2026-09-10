using CairoMakie: Figure
using Random
using Test

const _SCAN_EXAMPLE_PATH = joinpath(
    @__DIR__, "..", "..", "examples", "scan_rotated_planar.jl")
const _SCAN_CSV_HEADER =
    "distance,logical_state,construction,boundary_orientation,clock,p_x,p_z," *
    "shots,seed,logical_x_failures,logical_x_failure_rate," *
    "logical_x_standard_error,logical_z_failures,logical_z_failure_rate," *
    "logical_z_standard_error,any_logical_failures,any_logical_failure_rate," *
    "any_logical_standard_error,state_failures,state_failure_rate," *
    "state_standard_error"

@testset "Logical-failure scans and artifacts" begin
    scan_api = (
        :LogicalFailureScan,
        :scan_logical_failure,
        :write_logical_failure_csv,
        :plot_logical_failure_scan,
        :save_logical_failure_scan,
    )
    for name in scan_api
        @test isdefined(TopoNoise, name)
    end

    if all(name -> isdefined(TopoNoise, name), scan_api)
        @testset "axes are validated before simulation" begin
            invalid_distances = (Int[], [2], [3, 4], [3, 3], [5, 3])
            for distances in invalid_distances
                @test_throws ArgumentError scan_logical_failure(
                    MersenneTwister(1); distances, error_rates=[0.0], shots=1)
            end

            invalid_rates = (
                Float64[], [-0.01], [0.5], [Inf], [NaN],
                [0.0, 0.0], [0.1, 0.0],
            )
            for error_rates in invalid_rates
                @test_throws ArgumentError scan_logical_failure(
                    MersenneTwister(1);
                    distances=[3], error_rates, shots=1)
            end

            @test_throws ArgumentError scan_logical_failure(
                MersenneTwister(1);
                distances=[3], error_rates=[0.0], logical_state=:invalid, shots=1)
            @test_throws ArgumentError scan_logical_failure(
                MersenneTwister(1);
                distances=[3], error_rates=[0.0], construction=:invalid, shots=1)
            @test_throws ArgumentError scan_logical_failure(
                MersenneTwister(1);
                distances=[3], error_rates=[0.0], boundary_orientation=:invalid,
                shots=1)
            @test_throws ArgumentError scan_logical_failure(
                MersenneTwister(1);
                distances=[3], error_rates=[0.0], clock=:invalid, shots=1)
            @test_throws ArgumentError scan_logical_failure(
                MersenneTwister(1);
                distances=[3], error_rates=[0.0], shots=0)
            @test_throws ArgumentError scan_logical_failure(
                MersenneTwister(1);
                distances=[3], error_rates=[0.0], shots=1, batch_size=0)
        end

        @testset "Cartesian traversal is deterministic and preserves metadata" begin
            scan_rng = MersenneTwister(314)
            scan = scan_logical_failure(
                scan_rng;
                distances=[3, 5], error_rates=[0.0, 0.02],
                logical_state=:minus, construction=:as,
                boundary_orientation=:x_ew, clock=:plaquette,
                shots=2, batch_size=1, seed=314)

            @test scan isa LogicalFailureScan
            @test scan.distances == [3, 5]
            @test scan.error_rates == [0.0, 0.02]
            @test scan.logical_state === :minus
            @test scan.construction === :as
            @test scan.boundary_orientation === :x_ew
            @test scan.clock === :plaquette
            @test scan.shots == 2
            @test scan.batch_size == 1
            @test scan.seed == 314
            @test [(point.distance, point.p_x, point.p_z) for point in scan.points] == [
                (3, 0.0, 0.0), (3, 0.02, 0.02),
                (5, 0.0, 0.0), (5, 0.02, 0.02),
            ]
            @test all(point -> point.logical_state === :minus, scan.points)
            @test all(point -> point.construction === :as, scan.points)
            @test all(point -> point.boundary_orientation === :x_ew, scan.points)
            @test all(point -> point.clock === :plaquette, scan.points)
            @test all(point -> point.shots == 2, scan.points)
            @test all(point -> point.seed == 314, scan.points)

            manual_rng = MersenneTwister(314)
            manual_points = LogicalFailurePoint[]
            for d in (3, 5), p in (0.0, 0.02)
                code = RotatedPlanarCode(d; boundary_orientation=:x_ew)
                encoder = rotated_planar_encoder(
                    code; logical_state=:minus, construction=:as)
                push!(manual_points, estimate_logical_failure(
                    manual_rng, encoder, CircuitPauliNoise(p; clock=:plaquette);
                    shots=2, batch_size=1, seed=314))
            end
            @test [Tuple(getfield(point, field) for field in fieldnames(LogicalFailurePoint))
                   for point in scan.points] ==
                  [Tuple(getfield(point, field) for field in fieldnames(LogicalFailurePoint))
                   for point in manual_points]

            repeated = scan_logical_failure(
                MersenneTwister(314);
                distances=[3, 5], error_rates=[0.0, 0.02],
                logical_state=:minus, construction=:as,
                boundary_orientation=:x_ew, clock=:plaquette,
                shots=2, batch_size=1, seed=999)
            @test repeated.seed == 999
            for (original, relabeled) in zip(scan.points, repeated.points)
                for field in fieldnames(LogicalFailurePoint)
                    field === :seed ||
                        @test getfield(original, field) == getfield(relabeled, field)
                end
            end
        end

        @testset "zero noise, progress, and CSV contract" begin
            progress = IOBuffer()
            scan = scan_logical_failure(
                MersenneTwister(20);
                distances=[3], error_rates=[0.0],
                shots=2, batch_size=2, seed=nothing, progress_io=progress)
            @test length(scan.points) == 1
            @test scan.seed === nothing
            point = only(scan.points)
            @test point.logical_x_failures == 0
            @test point.logical_z_failures == 0
            @test point.any_logical_failures == 0
            @test point.state_failures == 0
            @test count(==('\n'), String(take!(progress))) == 1

            mktempdir() do directory
                empty_seed_path = joinpath(directory, "empty-seed.csv")
                @test write_logical_failure_csv(empty_seed_path, scan) ==
                    abspath(empty_seed_path)
                lines = readlines(empty_seed_path)
                @test lines[1] == _SCAN_CSV_HEADER
                @test lines[2] ==
                    "3,zero,bp,x_ns,gate_layer,0.0,0.0,2,,0,0.0,0.0," *
                    "0,0.0,0.0,0,0.0,0.0,0,0.0,0.0"

                seeded_scan = scan_logical_failure(
                    MersenneTwister(20);
                    distances=[3], error_rates=[0.0],
                    shots=2, batch_size=2, seed=20)
                seeded_path = joinpath(directory, "seeded.csv")
                write_logical_failure_csv(seeded_path, seeded_scan)
                @test split(readlines(seeded_path)[2], ',')[9] == "20"
            end
        end

        @testset "plot and four-file save smoke test" begin
            scan = scan_logical_failure(
                MersenneTwister(7);
                distances=[3], error_rates=[0.0],
                shots=2, batch_size=2, seed=7)
            @test plot_logical_failure_scan(scan) isa Figure

            mktempdir() do root
                output_dir = joinpath(root, "nested", "results")
                paths = save_logical_failure_scan(
                    scan, output_dir; basename="smoke")
                @test keys(paths) == (:csv, :svg, :pdf, :png)
                @test paths == (
                    csv=abspath(joinpath(output_dir, "smoke.csv")),
                    svg=abspath(joinpath(output_dir, "smoke.svg")),
                    pdf=abspath(joinpath(output_dir, "smoke.pdf")),
                    png=abspath(joinpath(output_dir, "smoke.png")),
                )
                for path in paths
                    @test isfile(path)
                    @test filesize(path) > 100
                end
                @test_throws ArgumentError save_logical_failure_scan(
                    scan, output_dir; basename="../escape")
            end
        end
    end
end

@testset "Rotated-planar scan CLI" begin
    @test isfile(_SCAN_EXAMPLE_PATH)

    if isfile(_SCAN_EXAMPLE_PATH)
        include(_SCAN_EXAMPLE_PATH)

        mktempdir() do directory
            output = IOBuffer()
            arguments = [
                "--distances", "3",
                "--error-rates", "0",
                "--shots", "2",
                "--batch-size", "2",
                "--seed", "41",
                "--output-dir", directory,
                "--basename", "cli-smoke",
            ]
            @test RotatedPlanarScanExample.main(
                arguments; io=output, error_io=output) == 0
            text = String(take!(output))
            for extension in ("csv", "svg", "pdf", "png")
                path = abspath(joinpath(directory, "cli-smoke.$extension"))
                @test isfile(path)
                @test filesize(path) > 100
                @test contains(text, path)
            end
        end

        output = IOBuffer()
        @test RotatedPlanarScanExample.main(
            ["--help"]; io=output, error_io=output) == 0
        @test contains(String(take!(output)), "Usage:")

        output = IOBuffer()
        @test RotatedPlanarScanExample.main(
            ["--unknown"]; io=output, error_io=output) == 1
        text = String(take!(output))
        @test contains(text, "error:")
        @test contains(text, "Usage:")
        @test contains(text, "unknown option")

        for arguments in (
                ["--p-max", "Inf"],
                ["--p-min", "-Inf"],
                ["--p-step", "1e-320"],
            )
            output = IOBuffer()
            outcome = try
                (
                    result=RotatedPlanarScanExample.main(
                        arguments; io=output, error_io=output),
                    error=nothing,
                )
            catch error
                (result=nothing, error)
            end
            @test outcome.error === nothing
            @test outcome.result == 1
            text = String(take!(output))
            @test contains(text, "error:")
            @test contains(text, "Usage:")
        end
    end
end
