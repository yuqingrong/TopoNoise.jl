@testset "Native product inputs pair every noisy H and CNOT" begin
    for d in (3, 5, 9, 11, 13, 15), orientation in (:x_ns, :x_ew)
        code = RotatedPlanarCode(d; boundary_orientation=orientation)
        as = rotated_planar_encoder(code; construction=:as)
        bp = rotated_planar_encoder(code; construction=:bp)
        rotation = [(column - 1) * d + d + 1 - row
                    for row in 1:d for column in 1:d]
        operations(encoder) = [op for layer in gate_layers(encoder)
                              for op in layer.operations]
        signature(encoder) = [(op.gate, op.qubits) for op in operations(encoder)]
        # Pair the whole executable schedule, including H faults; a check
        # of CNOTs alone misses the previously asymmetric preparation noise.
        expected = [(op.gate, op.gate === :CNOT ?
                     [rotation[op.qubits[2]], rotation[op.qubits[1]]] :
                     rotation[op.qubits]) for op in operations(bp)]
        @test signature(as) == expected
        @test hasproperty(as, :input_state)
        if hasproperty(as, :input_state)
            @test as.input_state === :plus
            @test bp.input_state === :zero
        end
        @test verify_encoder_tableau(as)
        @test verify_encoder_tableau(bp)
        @test sum(length(op.qubits) for op in operations(as)) ==
              sum(length(op.qubits) for op in operations(bp))
    end
end
