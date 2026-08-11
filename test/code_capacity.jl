# PyMatching imports Matplotlib lazily; keep its cache inside a writable test tempdir.
ENV["MPLCONFIGDIR"] = mktempdir()

using TopoNoise
using Test
using Random

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
