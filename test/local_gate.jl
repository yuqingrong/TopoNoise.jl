using LinearAlgebra: I

@testset "Unitary completion" begin
    coefficient = inv(sqrt(2))
    isometry = ComplexF64[
        coefficient 0
        coefficient * im 0
        0 1
        0 0
    ]

    unitary = unitary_completion(isometry)
    identity₄ = Matrix{ComplexF64}(I, 4, 4)
    @test size(unitary) == (4, 4)
    @test unitary[:, 1:2] == isometry
    @test adjoint(unitary) * unitary ≈ identity₄
    @test unitary * adjoint(unitary) ≈ identity₄

    square_unitary = ComplexF64[0 1; 1 0]
    @test unitary_completion(square_unitary) == square_unitary

    @test_throws DimensionMismatch unitary_completion(zeros(Float64, 2, 3))
    @test_throws ArgumentError unitary_completion(reshape([1, 0], 2, 1))
    @test_throws ArgumentError unitary_completion(reshape([1.0, 1.0], 2, 1))
    @test_throws ArgumentError unitary_completion(
        Matrix{Float64}(I, 2, 1); atol=-1)
    @test_throws ArgumentError unitary_completion(
        Matrix{Float64}(I, 2, 1); rtol=-1)
end

@testset "Toric-code local unitary gate" begin
    for scalar_type in (Float64, ComplexF64)
        tensor = toric_code_local_tensor(scalar_type)
        gate = toric_code_local_gate(scalar_type)
        gate_matrix = reshape(gate, 64, 64)
        identity₆₄ = Matrix{scalar_type}(I, 64, 64)

        @test size(gate) == ntuple(_ -> 2, 12)
        @test eltype(gate) == scalar_type
        @test adjoint(gate_matrix) * gate_matrix ≈ identity₆₄
        @test gate_matrix * adjoint(gate_matrix) ≈ identity₆₄

        # Julia index 1 represents qubit value 0. Fixing all four physical
        # input legs therefore gives the original rank-eight tensor exactly.
        @test gate[:, :, :, :, :, :, 1, 1, 1, 1, :, :] == tensor

        # Check the same map one computational-basis input at a time.
        for gamma in 1:2, delta in 1:2
            @test gate[:, :, :, :, :, :, 1, 1, 1, 1, gamma, delta] ==
                  tensor[:, :, :, :, :, :, gamma, delta]
        end
    end

    @test_throws ArgumentError toric_code_local_gate(Int)
end
