using Test
using TopoNoise

@testset "Diagonal slice transfer prepares all physical checkpoint checks" begin
    @test isdefined(TopoNoise, :Toric3DSliceEncoder)
    @test isdefined(TopoNoise, :toric3d_slice_encoder)
    if isdefined(TopoNoise, :toric3d_slice_encoder)
        for (steps, ly, lz) in ((1, 2, 2), (2, 2, 2), (3, 2, 3))
            code = TopoNoise.Toric3DSliceLattice(steps, ly, lz)
            encoder = TopoNoise.toric3d_slice_encoder(code)
            cells = ly * lz
            @test encoder.code === code
            @test encoder.checkpoints == [7 * cells * t for t in 1:steps]
            @test length(encoder.operations) == 7 * cells * steps
            @test count(op -> op.gate === :H, encoder.operations) == cells * steps
            @test count(op -> op.gate === :CNOT, encoder.operations) == 6 * cells * steps
            @test sum(length(op.qubits) for op in encoder.operations) == 13 * cells * steps
            @test TopoNoise.verify_toric3d_slice_encoder(encoder)
            touched = Set{Int}()
            for op in encoder.operations
                op.gate === :H && @test only(op.qubits) ∉ touched
                if op.gate === :CNOT
                    control = data_qubit_coordinate(code, op.qubits[1])
                    target = data_qubit_coordinate(code, op.qubits[2])
                    @test target[2] - control[2] in (0, 1)
                    control[2] == target[2] && @test control[1] == 1
                end
                union!(touched, op.qubits)
            end

            for t in 1:steps
                prefix = TopoNoise.Toric3DSliceLattice(t, ly, lz)
                operations = encoder.operations[1:encoder.checkpoints[t]]
                checks = ((:X, star_checks(prefix)),
                          (:Z, plaquette_checks(prefix)),
                          (:Z, TopoNoise.toric3d_slice_boundary_checks(prefix)),
                          (:Z, TopoNoise.toric3d_slice_checks(prefix, t)))
                @test all(all(data_qubit_coordinate(code, q)[2] <= t for q in op.qubits)
                          for op in operations)
                for (pauli, supports) in checks, support in supports
                    @test TopoNoise._toric3d_pulled_to_zero_input(
                        data_qubit_count(code), operations, pauli, support)
                end
            end

            shortened = TopoNoise.Toric3DSliceEncoder(code,
                encoder.operations[1:end-1], encoder.checkpoints)
            @test !TopoNoise.verify_toric3d_slice_encoder(shortened)
            invalid_checkpoints = TopoNoise.Toric3DSliceEncoder(code,
                encoder.operations, encoder.checkpoints .+ 1)
            @test !TopoNoise.verify_toric3d_slice_encoder(invalid_checkpoints)
            wrong = copy(encoder.operations)
            wrong[1] = TopoNoise.EncoderOperation(:CNOT, wrong[2].qubits)
            @test !TopoNoise.verify_toric3d_slice_encoder(
                TopoNoise.Toric3DSliceEncoder(code, wrong, encoder.checkpoints))
        end
        small = TopoNoise.toric3d_slice_encoder(TopoNoise.Toric3DSliceLattice(2, 2, 2))
        @test [op.gate for op in small.operations[1:7]] == [:H, :CNOT, :CNOT, :CNOT, :CNOT, :CNOT, :CNOT]
        @test [op.qubits for op in small.operations[1:7]] ==
              [[13], [13, 17], [13, 21], [7, 17], [2, 17], [11, 21], [1, 21]]
    end
end
