# PyMatching imports Matplotlib lazily; keep its cache inside a writable test tempdir.
ENV["MPLCONFIGDIR"] = mktempdir()

using TopoNoise
using Test
using Random

function synthetic_rotated_scan(
        error_rates::Vector{Float64}, small_rates::Vector{Float64},
        large_rates::Vector{Float64})
    length(error_rates) == length(small_rates) == length(large_rates) ||
        throw(ArgumentError("synthetic crossing curves must share a grid"))
    points = RotatedCodeCapacityScanPoint[]
    for (distance, rates) in ((3, small_rates), (5, large_rates))
        for (error_rate, failure_rate) in zip(error_rates, rates)
            push!(points, RotatedCodeCapacityScanPoint(
                distance, error_rate, 20, round(Int, 20 * failure_rate),
                failure_rate, 0.0, fill(failure_rate, 4), missing))
        end
    end
    return RotatedCodeCapacityScan([3, 5], error_rates, points)
end

@testset "rotated planar code-capacity circuit" begin
    @test_throws ArgumentError RotatedCodeCapacityModel(0)
    model = RotatedCodeCapacityModel(3)
    @test model.distance == 3
    @test length(model.data_qubits) == 9

    noisy = rotated_code_capacity_circuit(model, 0.1)
    @test occursin("X_ERROR(0.1)", string(noisy))
    @test all(occursin(string(qubit), string(noisy)) for qubit in model.data_qubits)

    clean = rotated_code_capacity_circuit(model, 0.0)
    dets, observables = clean.compile_detector_sampler().sample(
        shots=8, separate_observables=true)
    @test !any(TopoNoise.pyconvert(BitMatrix, dets))
    @test !any(TopoNoise.pyconvert(BitMatrix, observables))
    @test length(noisy.shortest_graphlike_error()) == 3
    @test TopoNoise._minimum_x_logical_support(model) == [1, 8, 15]
end

@testset "perfect Z-syndrome MWPM decoding" begin
    model = RotatedCodeCapacityModel(3)
    matching = TopoNoise._rotated_code_matching(model, 0.1)

    for qubit in model.data_qubits
        circuit = TopoNoise._deterministic_data_x_circuit(model, [qubit])
        syndrome_py, actual_py = circuit.compile_detector_sampler().sample(
            shots=1, separate_observables=true)
        syndrome = TopoNoise.pyconvert(BitMatrix, syndrome_py)
        actual = TopoNoise.pyconvert(BitMatrix, actual_py)
        predicted = TopoNoise.pyconvert(
            BitMatrix, TopoNoise._decode_observables(matching, syndrome_py))
        @test !any(predicted .!= actual)
    end

    support = TopoNoise._minimum_x_logical_support(model)
    circuit = TopoNoise._deterministic_data_x_circuit(model, support)
    syndrome_py, actual_py = circuit.compile_detector_sampler().sample(
        shots=1, separate_observables=true)
    syndrome = TopoNoise.pyconvert(BitMatrix, syndrome_py)
    actual = TopoNoise.pyconvert(BitMatrix, actual_py)
    @test !any(syndrome)
    @test only(vec(actual))
    @test_throws MethodError TopoNoise._decode_observables(
        matching, syndrome_py, actual_py)
end

@testset "single-point estimator" begin
    point = estimate_rotated_code_capacity(
        MersenneTwister(17), RotatedCodeCapacityModel(3), 0.0;
        shots=12, batches=3, seed=17)
    @test point.logical_x_failure_count == 0
    @test point.logical_x_failure_rate == 0.0
    @test length(point.logical_x_failure_batches) == 3
end

@testset "rotated finite-size scan" begin
    scan = scan_rotated_code_capacity(
        MersenneTwister(21), [3, 5], [0.0, 0.1];
        shots=12, batches=3, seed=21)
    @test length(scan.points) == 4
    @test all(point -> point.seed == 21, scan.points)
    @test all(point -> length(point.logical_x_failure_batches) == 3, scan.points)
    @test_throws ArgumentError scan_rotated_code_capacity(
        MersenneTwister(1), [3, 5], [0.1, 0.2]; shots=10, batches=3)
    @test_throws ArgumentError scan_rotated_code_capacity(
        MersenneTwister(1), [5, 3], [0.1, 0.2]; shots=12, batches=3)
    @test_throws ArgumentError scan_rotated_code_capacity(
        MersenneTwister(1), [3, 5], [0.2, 0.1]; shots=12, batches=3)
    @test_throws ArgumentError scan_rotated_code_capacity(
        MersenneTwister(1), [3, 5], [0.1, 0.5]; shots=12, batches=3)
end

@testset "rotated crossing statuses" begin
    @test only(estimate_rotated_code_crossings(
        MersenneTwister(2), synthetic_rotated_scan(
            [0.08, 0.10, 0.12], [0.2, 0.45, 0.7], [0.3, 0.45, 0.6]);
        bootstrap=20)).status == :ok
    @test only(estimate_rotated_code_crossings(
        MersenneTwister(3), synthetic_rotated_scan(
            [0.08, 0.10, 0.12], [0.2, 0.3, 0.4], [0.1, 0.2, 0.3]);
        bootstrap=20)).status == :unbracketed
    @test only(estimate_rotated_code_crossings(
        MersenneTwister(4), synthetic_rotated_scan(
            [0.08, 0.10, 0.12], [0.2, 0.5, 0.8], [0.2, 0.5, 0.8]);
        bootstrap=20)).status == :unstable

    unequal_batches = synthetic_rotated_scan(
        [0.08, 0.10, 0.12], [0.2, 0.45, 0.7], [0.3, 0.45, 0.6])
    pop!(first(unequal_batches.points).logical_x_failure_batches)
    @test_throws ArgumentError estimate_rotated_code_crossings(
        MersenneTwister(5), unequal_batches; bootstrap=20)

    complete_scan = synthetic_rotated_scan(
        [0.08, 0.10, 0.12], [0.2, 0.45, 0.7], [0.3, 0.45, 0.6])
    incomplete_scan = RotatedCodeCapacityScan(
        complete_scan.distances, complete_scan.error_rates,
        complete_scan.points[1:(end - 1)])
    @test_throws ArgumentError estimate_rotated_code_crossings(
        MersenneTwister(6), incomplete_scan; bootstrap=20)
end
