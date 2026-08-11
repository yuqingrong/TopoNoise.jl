using TopoNoise
using Test

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
