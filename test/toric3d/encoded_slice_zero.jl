using Test
using TopoNoise

function toric3d_zero_test_rank(matrix)
    work = BitMatrix(matrix)
    pivot = 0
    for column in axes(work, 2)
        source = findfirst(row -> work[row, column], pivot + 1:size(work, 1))
        isnothing(source) && continue
        source += pivot
        pivot += 1
        work[pivot, :], work[source, :] = copy(work[source, :]), copy(work[pivot, :])
        for row in pivot + 1:size(work, 1)
            work[row, column] && (work[row, :] .⊻= work[pivot, :])
        end
    end
    return pivot
end

@testset "Encoded transfer prepares logical zero at every checkpoint" begin
    for dims in ((1, 2, 2), (2, 2, 3), (3, 3, 3))
        code = Toric3DEncodedSlice(dims...)
        encoder = toric3d_encoded_encoder(code)
        @test TopoNoise.verify_toric3d_encoded_zero_encoder(code)
        @test TopoNoise.verify_toric3d_encoded_zero_encoder(code, encoder)

        short = TopoNoise.Toric3DSliceEncoder(code.lattice,
            encoder.operations[1:end-1], encoder.checkpoints)
        @test !TopoNoise.verify_toric3d_encoded_zero_encoder(code, short)

        shifted = TopoNoise.Toric3DSliceEncoder(code.lattice,
            encoder.operations, encoder.checkpoints .+ 1)
        @test !TopoNoise.verify_toric3d_encoded_zero_encoder(code, shifted)

        changed = copy(encoder.operations)
        changed[1] = TopoNoise.EncoderOperation(:CNOT, changed[2].qubits)
        corrupted = TopoNoise.Toric3DSliceEncoder(code.lattice,
            changed, encoder.checkpoints)
        @test !TopoNoise.verify_toric3d_encoded_zero_encoder(code, corrupted)

        m = dims[2] * dims[3]
        seed = toric3d_encoded_seed_checks(code)
        seed_x = TopoNoise._check_matrix(seed.x, 3m; sparse=false)
        seed_z = TopoNoise._check_matrix(seed.z, 3m; sparse=false)
        logical_z = TopoNoise._check_matrix(toric3d_encoded_logical_z(code), 3m; sparse=false)
        @test toric3d_zero_test_rank(seed_x) == 2m - 1
        @test toric3d_zero_test_rank(seed_z) == m - 1
        @test toric3d_zero_test_rank(vcat(seed_z, logical_z)) == m + 1
    end

    periodic_two = Toric3DEncodedSlice(1, 2, 2)
    n = data_qubit_count(periodic_two)
    logical_x = TopoNoise._check_matrix(toric3d_encoded_logical_x(periodic_two), n;
        sparse=false)
    logical_z = TopoNoise._check_matrix(toric3d_encoded_logical_z(periodic_two), n;
        sparse=false)
    @test Int.(logical_x) * transpose(Int.(logical_z)) .% 2 == [1 0; 0 1]
end
