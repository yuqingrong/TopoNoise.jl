using Test
using TopoNoise

@testset "Encoded diagonal prefixes preserve full-code ordering" begin
    code = Toric3DEncodedSlice(3, 2, 3)
    for t in 1:3
        prefix = toric3d_encoded_prefix(code, t)
        @test prefix isa Toric3DEncodedPrefix
        @test (prefix.code.lattice.steps, prefix.code.lattice.Ly, prefix.code.lattice.Lz) ==
              (t, 2, 3)
        @test length(prefix.qubits) == data_qubit_count(prefix.code)
        @test length(unique(prefix.qubits)) == length(prefix.qubits)
        @test all(data_qubit_coordinate(prefix.code, q) ==
                  data_qubit_coordinate(code, prefix.qubits[q])
                  for q in eachindex(prefix.qubits))
        @test all(data_qubit_coordinate(code, q)[2] <= t for q in prefix.qubits)
    end
    @test_throws ArgumentError toric3d_encoded_prefix(code, 0)
    @test_throws ArgumentError toric3d_encoded_prefix(code, 4)
end

@testset "Z feedback windows are syndrome complete and constant width" begin
    for dims in ((2, 2, 2), (2, 2, 3), (3, 3, 2))
        code = Toric3DEncodedSlice(dims...)
        T, Ly, Lz = dims
        cells = Ly * Lz

        seed = toric3d_z_feedback_window(code, 0)
        @test seed isa Toric3DZFeedbackWindow
        @test seed.r == 0
        @test size(seed.matrix) == (2cells, 3cells)
        @test seed.qubits == TopoNoise.toric3d_slice_qubits(code.lattice, 0)
        @test Set(data_qubit_coordinate(code, q)[2] for check in seed.checks for q in check) ==
              Set([0])
        @test all(==(2), vec(sum(seed.matrix; dims=1)))

        for r in 1:T
            window = toric3d_z_feedback_window(code, r)
            @test window.r == r
            @test size(window.matrix) == (3cells, 6cells)
            @test window.qubits == vcat(
                TopoNoise.toric3d_slice_qubits(code.lattice, r - 1),
                TopoNoise.toric3d_slice_qubits(code.lattice, r))
            @test length(unique(window.qubits)) == 6cells
            @test all(==(2), vec(sum(window.matrix; dims=1)))
            measured_slices = Set(data_qubit_coordinate(code, q)[2]
                                  for check in window.checks for q in check)
            @test measured_slices == Set(r == 1 ? (0:1) : (r - 2:r))

            frame = falses(data_qubit_count(code))
            @test toric3d_z_window_syndrome(window, frame) == falses(3cells)
            frame[window.qubits[1]] = true
            @test toric3d_z_window_syndrome(window, frame) == window.matrix[:, 1]
            @test_throws DimensionMismatch toric3d_z_window_syndrome(
                window, falses(data_qubit_count(code) - 1))
        end

        @test_throws ArgumentError toric3d_z_feedback_window(code, -1)
        @test_throws ArgumentError toric3d_z_feedback_window(code, T + 1)
    end
end

@testset "Length-two transverse directions retain parallel physical edges" begin
    code = Toric3DEncodedSlice(2, 2, 2)
    window = toric3d_z_feedback_window(code, 1)
    @test length(window.qubits) == 24
    @test length(unique(window.qubits)) == 24
    @test length(unique(data_qubit_coordinate(code, q) for q in window.qubits)) == 24
    periodic_pairs = (
        (data_qubit_index(code, 1, 0, 0, 0), data_qubit_index(code, 1, 0, 1, 0)),
        (data_qubit_index(code, 2, 0, 0, 0), data_qubit_index(code, 2, 0, 0, 1)),
    )
    for (first_edge, second_edge) in periodic_pairs
        first_column = findfirst(==(first_edge), window.qubits)
        second_column = findfirst(==(second_edge), window.qubits)
        @test !isnothing(first_column)
        @test !isnothing(second_column)
        @test first_edge != second_edge
        @test window.matrix[:, first_column] != window.matrix[:, second_column]
    end
end
