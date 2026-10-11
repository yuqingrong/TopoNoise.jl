using Test
using TopoNoise

function toric3d_test_pullback(n, operations, pauli, support)
    x, z = falses(n), falses(n)
    (pauli === :X ? x : z)[support] .= true
    for operation in Iterators.reverse(operations)
        if operation.gate === :H
            q = only(operation.qubits)
            x[q], z[q] = z[q], x[q]
        else
            control, target = operation.qubits
            x[target] ⊻= x[control]
            z[control] ⊻= z[target]
        end
    end
    return x, z
end

@testset "3D sequential encoder" begin
    code = ToricCode3D(2, 2, 2)
    encoder = toric3d_encoder(code)
    @test length(encoder.checkpoints) == 2
    @test length(encoder.operations) == 6 * 2 * 2 * 2
    @test encoder.checkpoints == [24, 48]
    @test verify_toric3d_encoder(encoder)
    @test all(op.gate === :H || op.gate === :CNOT for op in encoder.operations)
    @test all(op.gate === :H for op in encoder.operations[1:6:48])
    touched = Set{Int}()
    for operation in encoder.operations
        operation.gate === :H && @test only(operation.qubits) ∉ touched
        union!(touched, operation.qubits)
    end
    for j in 1:2
        prefix = encoder.prefixes[j]
        mapping = encoder.prefix_qubits[j]
        checkpoint = encoder.checkpoints[j]
        used = Set{Int}()
        for operation in encoder.operations[1:checkpoint]
            union!(used, operation.qubits)
        end
        @test used == Set(mapping)
        @test all(all(q in mapping for q in operation.qubits)
                  for operation in encoder.operations[1:checkpoint])
        @test all(data_qubit_index(code, :x, x, y, z) ∉ used
                  for x in j+1:code.Lx for y in 0:code.Ly-1 for z in 0:code.Lz-1)
        for (pauli, checks) in ((:X, star_checks(prefix)),
                                (:Z, plaquette_checks(prefix)),
                                (:Z, [logical_z_support(prefix)]))
            for support in checks
                x, _ = toric3d_test_pullback(data_qubit_count(code),
                    encoder.operations[1:checkpoint], pauli, mapping[support])
                @test !any(x)
            end
        end
    end
    @test encoder.operations[1].qubits ==
          [data_qubit_index(code, :x, 1, 0, 0)]
    @test [op.qubits[2] for op in encoder.operations[2:6]] == [
        data_qubit_index(code, :x, 0, 0, 0),
        data_qubit_index(code, :y, 1, 0, 0),
        data_qubit_index(code, :y, 1, 1, 0),
        data_qubit_index(code, :z, 1, 0, 0),
        data_qubit_index(code, :z, 1, 0, 1),
    ]
    @test !verify_toric3d_encoder(Toric3DEncoder(code,
        encoder.operations[1:end-1], encoder.checkpoints,
        encoder.prefixes, encoder.prefix_qubits))
end
