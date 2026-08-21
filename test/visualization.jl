using CairoMakie: CairoMakie, Figure
using Random
using Test
using TopoNoise
import Yao

function rendered_text(fig)
    labels = String[]
    for block in fig.content
        hasproperty(block, :scene) || continue
        for plot in block.scene.plots
            hasproperty(plot, :text) || continue
            append!(labels, string.(plot.text[]))
        end
    end
    return labels
end

@testset "Rotated code-capacity figure" begin
    @test isdefined(TopoNoise, :plot_rotated_code_capacity)

    scan = RotatedCodeCapacityScan(
            [3, 5, 7], [0.05, 0.10], [
                RotatedCodeCapacityScanPoint(3, 0.05, 12, 1, 1 / 12,
                    sqrt((1 / 12) * (11 / 12) / 12), [0.0, 0.0, 0.25], 31),
                RotatedCodeCapacityScanPoint(3, 0.10, 12, 3, 0.25,
                    sqrt(0.25 * 0.75 / 12), [0.0, 0.25, 0.5], 31),
                RotatedCodeCapacityScanPoint(5, 0.05, 12, 2, 1 / 6,
                    sqrt((1 / 6) * (5 / 6) / 12), [0.0, 0.25, 0.25], 31),
                RotatedCodeCapacityScanPoint(5, 0.10, 12, 3, 0.25,
                    sqrt(0.25 * 0.75 / 12), [0.25, 0.25, 0.25], 31),
                RotatedCodeCapacityScanPoint(7, 0.05, 12, 3, 0.25,
                    sqrt(0.25 * 0.75 / 12), [0.25, 0.25, 0.25], 31),
                RotatedCodeCapacityScanPoint(7, 0.10, 12, 6, 0.5,
                    sqrt(0.5 * 0.5 / 12), [0.25, 0.5, 0.75], 31),
            ])
        crossings = [
            CriticalCrossing(3, 5, 0.075, 0.060, 0.090, 0.95, 1.0, :ok),
            CriticalCrossing(5, 7, missing, missing, missing, 0.95, 0.0,
                :unbracketed),
        ]
    figure = plot_rotated_code_capacity(scan, crossings)
    labels = join(rendered_text(figure), "\n")

    @test figure isa Figure
    @test occursin("p_c", labels)
    @test occursin("unbracketed", labels)
end

@testset "Isometric planar capacity scaling figure" begin
    @test isdefined(TopoNoise, :fit_isometric_planar_scaling)
    @test isdefined(TopoNoise, :plot_isometric_planar_scaling_diagnostics)

    if isdefined(TopoNoise, :fit_isometric_planar_scaling) &&
       isdefined(TopoNoise, :plot_isometric_planar_scaling_diagnostics)
        scan = IsometricPlanarCapacityScan(
            [5, 7, 9], [0.08, 0.10, 0.12], [
                IsometricPlanarCapacityPoint(
                    distance, rate, 400, round(Int, 400 * (0.25 + 4 * (rate - 0.10))),
                    0.25 + 4 * (rate - 0.10), 0.01,
                    fill(0.25 + 4 * (rate - 0.10), 4), false, 9)
                for distance in [5, 7, 9] for rate in [0.08, 0.10, 0.12]
            ])
        fit = fit_isometric_planar_scaling(
            MersenneTwister(73), scan; bootstrap=12)
        crossings = [
            CriticalCrossing(5, 7, 0.095, 0.09, 0.10, 0.95, 1.0, :ok),
            CriticalCrossing(7, 9, 0.105, 0.10, 0.11, 0.95, 1.0, :ok),
        ]
        figure = plot_isometric_planar_capacity(
            scan, crossings; scaling_fit=fit)
        @test figure isa Figure
        capacity_titles = [string(block.title[]) for block in figure.content
                           if hasproperty(block, :title)]
        @test !("Adjacent-distance crossing summary" in capacity_titles)
        capacity_axis = only(block for block in figure.content
                             if hasproperty(block, :title) &&
                                string(block.title[]) ==
                                "Isometric planar-code logical failure")
        @test capacity_axis.limits[][2] == (-0.01, 0.51)
        @test capacity_axis.yticks[] == 0.0:0.1:0.5
        @test !capacity_axis.xgridvisible[]
        @test !capacity_axis.ygridvisible[]
        line_plots = [plot for plot in capacity_axis.scene.plots
                      if occursin("Lines", string(typeof(plot)))]
        @test length(line_plots) == length(scan.distances) + 1
        @test !any(plot -> occursin("VLines", string(typeof(plot))),
                   capacity_axis.scene.plots)
        scaling_line_points = last(line_plots)[1][]
        @test length(scaling_line_points) == 2
        @test all(point -> all(isfinite, point), scaling_line_points)
        @test first(scaling_line_points)[2] == -0.01
        @test last(scaling_line_points)[2] == 0.51

        diagnostics = diagnose_isometric_planar_scaling(
            MersenneTwister(75), scan; bootstrap=12,
            loss_grid_size=(11, 9), sensitivity=true)
        diagnostic_figure = plot_isometric_planar_scaling_diagnostics(
            scan, diagnostics)
        @test diagnostic_figure isa Figure
        diagnostic_titles = [string(block.title[]) for block in diagnostic_figure.content
                             if hasproperty(block, :title)]
        @test "Profiled loss surface" in diagnostic_titles
        @test "Bootstrap parameter distribution" in diagnostic_titles
        @test "Batch variation check" in diagnostic_titles
        @test "Sensitivity: p_c" in diagnostic_titles
        @test "Sensitivity: ν" in diagnostic_titles
        diagnostic_labels = join(rendered_text(diagnostic_figure), "\n")
        @test occursin("degenerate bootstrap", diagnostic_labels)
        loss_axis = only(block for block in diagnostic_figure.content
                         if hasproperty(block, :title) &&
                            string(block.title[]) == "Profiled loss surface")
        @test length(loss_axis.scene.plots) >= 4

        short_fit = fit_isometric_planar_scaling(
            MersenneTwister(74), IsometricPlanarCapacityScan(
                [5, 7], scan.error_rates,
                filter(point -> point.distance in (5, 7), scan.points)); bootstrap=12)
        short_figure = plot_isometric_planar_capacity(
            IsometricPlanarCapacityScan(
                [5, 7], scan.error_rates,
                filter(point -> point.distance in (5, 7), scan.points)),
            CriticalCrossing[]; scaling_fit=short_fit)
        @test any(occursin("fit unavailable", label) for label in rendered_text(short_figure))

        short_diagnostics = diagnose_isometric_planar_scaling(
            MersenneTwister(76), IsometricPlanarCapacityScan(
                [5, 7], scan.error_rates,
                filter(point -> point.distance in (5, 7), scan.points));
            bootstrap=12, loss_grid_size=(11, 9), sensitivity=true)
        short_diagnostic_figure = plot_isometric_planar_scaling_diagnostics(
            IsometricPlanarCapacityScan(
                [5, 7], scan.error_rates,
                filter(point -> point.distance in (5, 7), scan.points)),
            short_diagnostics)
        @test any(occursin("diagnostics unavailable", label)
                  for label in rendered_text(short_diagnostic_figure))
    end
end

@testset "Trajectory scan figure" begin
    @test isdefined(TopoNoise, :plot_trajectory_scan)

    if isdefined(TopoNoise, :plot_trajectory_scan)
        scan = scan_trajectories(
            MersenneTwister(301), [2, 3], [0.0, 0.5, 1.0];
            shots=12, batches=3)
        crossings = estimate_crossings(
            MersenneTwister(302), scan; bootstrap=12)
        binder_crossings = estimate_binder_crossings(
            MersenneTwister(303), scan; bootstrap=12)
        figure = plot_trajectory_scan(
            scan, crossings; binder_crossings=binder_crossings)
        @test figure isa Figure

        titles = [string(block.title[]) for block in figure.content
                  if hasproperty(block, :title)]
        @test "Injected errors and mismatches" in titles
        @test "Plaquette frustration" in titles
        @test "Largest occupied cluster" in titles
        @test "Horizontal spanning probability" in titles
        @test "Marginalized |M|" in titles
        @test "Binder cumulant U₄" in titles
        @test !("Logical failure (N-S, MWPM)" in titles)

        raw_only_scan = TrajectoryScan(
            scan.sizes, scan.error_rates, scan.points)
        raw_only_figure = plot_trajectory_scan(raw_only_scan, crossings)
        @test raw_only_figure isa Figure
        raw_only_titles = [
            string(block.title[]) for block in raw_only_figure.content
            if hasproperty(block, :title)]
        @test !("Marginalized |M|" in raw_only_titles)
        @test !("Binder cumulant U₄" in raw_only_titles)

        negative_spin_points = [
            MarginalSpinScanPoint(
                point.size, point.error_rate, point.shots,
                point.absolute_magnetization_mean,
                point.absolute_magnetization_se,
                point.second_moment_mean, point.second_moment_se,
                point.fourth_moment_mean, point.fourth_moment_se,
                index == 1 ? -0.4 : point.binder_cumulant,
                missing,
                point.batch_counts, point.second_moment_batches,
                point.fourth_moment_batches)
            for (index, point) in enumerate(scan.marginal_spin_points)]
        negative_scan = TrajectoryScan(
            scan.sizes, scan.error_rates, scan.points, negative_spin_points)
        negative_figure = plot_trajectory_scan(negative_scan, crossings)
        negative_binder_axis = only(
            block for block in negative_figure.content
            if hasproperty(block, :title) &&
               string(block.title[]) == "Binder cumulant U₄")
        @test negative_binder_axis.limits[][2][1] < -0.4
        @test !any(occursin("Errorbars", string(typeof(plot)))
                   for plot in negative_binder_axis.scene.plots)

        mktempdir() do directory
            for extension in ("svg", "pdf", "png")
                path = joinpath(directory, "trajectory_scan.$extension")
                CairoMakie.save(path, figure)
                @test isfile(path)
                @test filesize(path) > 100
            end
            @test startswith(
                read(joinpath(directory, "trajectory_scan.svg"), String),
                "<?xml")
            @test startswith(
                read(joinpath(directory, "trajectory_scan.pdf"), String),
                "%PDF")
        end
    end
end

@testset "Native Yao circuit figures" begin
    @test isdefined(TopoNoise, :yao_circuit)

    if isdefined(TopoNoise, :yao_circuit)
        mktempdir() do directory
            for (rows, cols) in ((1, 1), (2, 2))
                block = yao_circuit(toric_code_peps(rows, cols))
                @test Yao.plot(block) !== nothing

                for format in (:svg, :pdf)
                    path = joinpath(directory, "yao-$(rows)x$(cols).$(format)")
                    Yao.vizcircuit(block; format=format, filename=path)
                    @test isfile(path)
                    @test filesize(path) > 100
                end
            end
            @test startswith(
                read(joinpath(directory, "yao-1x1.svg"), String), "<?xml")
            @test startswith(
                read(joinpath(directory, "yao-2x2.pdf"), String), "%PDF")
        end
    end
end

@testset "Toric-code vector figures" begin
    @test isdefined(TopoNoise, :plot_peps_graph)
    @test isdefined(TopoNoise, :plot_sequential_circuit)

    if isdefined(TopoNoise, :plot_peps_graph)
        peps = toric_code_peps(2, 2)
        circuit = sequential_circuit_graph(peps)
        peps_figure = plot_peps_graph(peps; show_index_labels=true)
        circuit_figure = plot_sequential_circuit(circuit)
        expanded_figure = plot_sequential_circuit(
            circuit; show_layer_labels=false, expand_physical_buses=true)

        @test peps_figure isa Figure
        @test circuit_figure isa Figure
        @test expanded_figure isa Figure

        peps_labels = rendered_text(peps_figure)
        @test all("T[$row,$col]" in peps_labels for row in 1:2 for col in 1:2)
        @test all(any(startswith("$direction:"), peps_labels)
                  for direction in ("E", "N", "W", "S"))

        circuit_labels = rendered_text(circuit_figure)
        @test all("U[$row,$col]" in circuit_labels for row in 1:2 for col in 1:2)
        @test any(contains("W→E"), circuit_labels)
        @test any(contains("S→N"), circuit_labels)
        @test "|0⟩" in circuit_labels
        @test "|+⟩" in circuit_labels
        @test "⟨+|" in circuit_labels
        @test all("layer $layer" in circuit_labels for layer in 1:3)

        mktempdir() do directory
            outputs = (
                joinpath(directory, "peps.svg") => peps_figure,
                joinpath(directory, "peps.pdf") => peps_figure,
                joinpath(directory, "circuit.svg") => circuit_figure,
                joinpath(directory, "circuit.pdf") => circuit_figure,
            )
            for (path, figure) in outputs
                CairoMakie.save(path, figure)
                @test isfile(path)
                @test filesize(path) > 100
            end
            @test startswith(read(joinpath(directory, "peps.svg"), String), "<?xml")
            @test startswith(read(joinpath(directory, "circuit.pdf"), String), "%PDF")
        end
    end
end
