using TopoNoise
using ITensors
using Test

function ordered_physical_indices(peps)
    rows, cols = size(peps)
    return [index for row in 1:rows for col in 1:cols
            for index in physicalinds(peps, row, col)]
end

function contract_state(peps)
    ITensors.disable_warn_order()
    try
        return reduce(*, vec(peps.tensors))
    finally
        ITensors.reset_warn_order()
    end
end

function position(cols, row, col, direction)
    return ((row - 1) * cols + col - 1) * 4 + direction
end

function is_allowed_configuration(configuration, rows, cols)
    bit(row, col, direction) =
        configuration[position(cols, row, col, direction)] - 1

    for row in 1:rows, col in 1:cols
        iseven(sum(bit(row, col, direction) for direction in 1:4)) || return false
    end
    for row in 1:rows, col in 1:(cols - 1)
        bit(row, col, 1) == bit(row, col + 1, 3) || return false
    end
    for row in 1:(rows - 1), col in 1:cols
        bit(row, col, 4) == bit(row + 1, col, 2) || return false
    end
    return true
end

function projected_amplitude(peps, bits)
    rows, cols = size(peps)
    size(bits) == (rows, cols, 4) || throw(DimensionMismatch("invalid bit array"))
    projected = Matrix{ITensor}(undef, rows, cols)
    for row in 1:rows, col in 1:cols
        tensor = peps[row, col]
        for (direction, physical) in enumerate(physicalinds(peps, row, col))
            tensor *= onehot(physical => bits[row, col, direction] + 1)
        end
        projected[row, col] = tensor
    end
    return scalar(reduce(*, vec(projected)))
end

function double_layer_norm2(peps)
    rows, cols = size(peps)
    double_layer = Matrix{ITensor}(undef, rows, cols)
    for row in 1:rows, col in 1:cols
        ket = peps[row, col]
        physical = physicalinds(peps, row, col)
        virtual = Tuple(index for index in inds(ket) if index ∉ physical)
        bra = dag(prime(ket; inds=virtual))
        double_layer[row, col] = bra * ket
    end
    return scalar(reduce(*, vec(double_layer)))
end

@testset "Local toric-code tensor" begin
    tensor = toric_code_local_tensor()
    @test size(tensor) == ntuple(_ -> 2, 8)
    @test count(!iszero, tensor) == 8
    @test_throws ArgumentError toric_code_local_tensor(Int)

    for index in CartesianIndices(tensor)
        physical = Tuple(index[direction] - 1 for direction in 1:4)
        virtual = Tuple(index[direction] - 1 for direction in 5:8)
        allowed = physical == virtual && iseven(sum(virtual))
        expected = allowed ? inv(sqrt(2)) : 0.0
        @test tensor[index] ≈ expected
    end

    # Contract each physical leg with the uniform covector (1, 1). The copy
    # constraints identify those physical values with their corresponding
    # virtual values, leaving the even-parity tensor on the virtual legs.
    virtual_tensor = dropdims(sum(tensor; dims=(1, 2, 3, 4)); dims=(1, 2, 3, 4))
    W = reshape(virtual_tensor, 4, 4)
    expected_W = [
        1 0 0 1
        0 1 1 0
        0 1 1 0
        1 0 0 1
    ] / sqrt(2)
    @test W ≈ expected_W

    # Stored order: (pE, pN, pW, pS, vE, vN, vW, vS).
    # Rows are (γ, δ, i, j, k, l) and columns are (α, β).
    isometry = reshape(permutedims(tensor, (7, 8, 1, 2, 3, 4, 5, 6)), 64, 4)
    identity_αβ = [
        1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
    ]
    @test adjoint(isometry) * isometry ≈ identity_αβ

    @test eltype(toric_code_local_tensor(ComplexF64)) == ComplexF64
end

@testset "Finite PEPS structure" begin
    @test_throws ArgumentError toric_code_peps(0, 2)
    @test_throws ArgumentError toric_code_peps(2, 0)
    @test_throws ArgumentError toric_code_peps(2, 2; scalar_type=Int)

    peps = toric_code_peps(3, 3)
    @test size(peps) == (3, 3)
    @test size(peps, 1) == 3
    @test length(inds(peps[1, 1])) == 6
    @test length(inds(peps[1, 2])) == 7
    @test length(inds(peps[2, 2])) == 8
    @test length(unique(ordered_physical_indices(peps))) == 36

    internal_bonds = Index[]
    for row in 1:3, col in 1:2
        shared = collect(commoninds(peps[row, col], peps[row, col + 1]))
        @test length(shared) == 1
        append!(internal_bonds, shared)
    end
    for row in 1:2, col in 1:3
        shared = collect(commoninds(peps[row, col], peps[row + 1, col]))
        @test length(shared) == 1
        append!(internal_bonds, shared)
    end
    @test length(internal_bonds) == 12
    @test length(unique(internal_bonds)) == 12
    @test isempty(commoninds(peps[1, 1], peps[3, 3]))

    complex_peps = toric_code_peps(1, 1; scalar_type=ComplexF64)
    @test eltype(complex_peps[1, 1]) == ComplexF64
end

@testset "Complete small-patch wavefunctions" begin
    for (rows, cols) in ((1, 1), (1, 2), (2, 2))
        peps = toric_code_peps(rows, cols)
        state = contract_state(peps)
        physical = ordered_physical_indices(peps)
        amplitudes = Array(state, physical...)

        vertices = rows * cols
        internal_edges = rows * (cols - 1) + (rows - 1) * cols
        boundary_edges = 2 * rows + 2 * cols
        exponent = internal_edges + boundary_edges - vertices
        expected_amplitude = exp2(-exponent / 2)

        @test count(!iszero, amplitudes) == 2^exponent
        @test all(CartesianIndices(amplitudes)) do configuration
            expected = is_allowed_configuration(configuration, rows, cols) ?
                       expected_amplitude : 0.0
            isapprox(amplitudes[configuration], expected; atol=1e-14, rtol=1e-14)
        end
        @test isapprox(sum(abs2, amplitudes), 1; atol=1e-13, rtol=1e-13)
    end
end

@testset "3×3 patch" begin
    peps = toric_code_peps(3, 3)
    expected_amplitude = exp2(-15 / 2)

    all_zero = zeros(Int, 3, 3, 4)
    @test projected_amplitude(peps, all_zero) ≈ expected_amplitude

    plaquette_loop = copy(all_zero)
    plaquette_loop[1, 1, 1] = plaquette_loop[1, 2, 3] = 1
    plaquette_loop[1, 2, 4] = plaquette_loop[2, 2, 2] = 1
    plaquette_loop[2, 1, 1] = plaquette_loop[2, 2, 3] = 1
    plaquette_loop[1, 1, 4] = plaquette_loop[2, 1, 2] = 1
    @test projected_amplitude(peps, plaquette_loop) ≈ expected_amplitude

    boundary_string = copy(all_zero)
    boundary_string[1, 1, 2] = 1
    boundary_string[1, 1, 4] = boundary_string[2, 1, 2] = 1
    boundary_string[2, 1, 4] = boundary_string[3, 1, 2] = 1
    boundary_string[3, 1, 4] = 1
    @test projected_amplitude(peps, boundary_string) ≈ expected_amplitude

    odd_vertex = copy(all_zero)
    odd_vertex[1, 1, 2] = 1
    @test iszero(projected_amplitude(peps, odd_vertex))

    mismatched_copies = copy(all_zero)
    mismatched_copies[1, 1, 2] = 1
    mismatched_copies[1, 1, 1] = 1
    @test iszero(projected_amplitude(peps, mismatched_copies))

    @test isapprox(double_layer_norm2(peps), 1; atol=1e-12, rtol=1e-12)
end
