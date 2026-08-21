using Random
using Statistics: mean
import Yao

function test_virtual_errors(
        rows, cols;
        horizontal=falses(rows, cols - 1),
        vertical=falses(rows - 1, cols),
        west=falses(rows), east=falses(rows),
        south=falses(cols), north=falses(cols))
    return VirtualBondErrors(
        BitMatrix(horizontal), BitMatrix(vertical),
        BitVector(west), BitVector(east),
        BitVector(south), BitVector(north))
end

function measurement_key(measurements)
    rows, cols, _ = size(measurements)
    return Tuple(measurements[row, col, direction]
                 for row in 1:rows for col in 1:cols for direction in 1:4)
end

function exact_classical_distribution(model, errors)
    rows, cols = model.rows, model.cols
    sites = [site for layer in model.layers for site in layer]
    latent_bits = rows + cols + length(sites)
    probability = exp2(-latent_bits)
    distribution = Dict{Tuple,Float64}()
    for configuration in 0:(2^latent_bits - 1)
        bit(index) = Bool((configuration >> (index - 1)) & 1)
        horizontal_carriers = [bit(row) for row in 1:rows]
        vertical_carriers = [bit(rows + col) for col in 1:cols]
        measurements = Array{Bool}(undef, rows, cols, 4)
        for (site_number, (row, col)) in enumerate(sites)
            gamma = horizontal_carriers[row]
            delta = vertical_carriers[col]
            alpha = bit(rows + cols + site_number)
            beta = xor(alpha, gamma, delta)
            measurements[row, col, :] .= (alpha, beta, gamma, delta)
            horizontal_carriers[row] = col < cols ?
                xor(alpha, errors.horizontal_internal[row, col]) : alpha
            vertical_carriers[col] = row > 1 ?
                xor(beta, errors.vertical_internal[row - 1, col]) : beta
        end
        key = measurement_key(measurements)
        distribution[key] = get(distribution, key, 0.0) + probability
    end
    return distribution
end

function exact_yao_distribution(model, errors)
    rows, cols = model.rows, model.cols
    physical_count = 4 * rows * cols
    first_horizontal = physical_count + 1
    first_vertical = first_horizontal + rows
    wire_count = physical_count + rows + cols
    register = Yao.zero_state(wire_count)
    Yao.apply!(register, Yao.repeat(
        wire_count, Yao.H, first_horizontal:wire_count))
    for row in 1:rows
        errors.west_boundary[row] && Yao.apply!(register,
            Yao.repeat(wire_count, Yao.X, (first_horizontal + row - 1,)))
    end
    for col in 1:cols
        errors.south_boundary[col] && Yao.apply!(register,
            Yao.repeat(wire_count, Yao.X, (first_vertical + col - 1,)))
    end

    local_matrix = reshape(toric_code_local_gate(ComplexF64), 64, 64)
    for layer in model.layers, (row, col) in layer
        physical_wires = TopoNoise._physical_wires(cols, row, col)
        horizontal_wire = first_horizontal + row - 1
        vertical_wire = first_vertical + col - 1
        gate = Yao.subroutine(
            wire_count,
            TopoNoise._yao_local_gate(local_matrix, (row, col)),
            (physical_wires..., horizontal_wire, vertical_wire))
        Yao.apply!(register, gate)
        col < cols && errors.horizontal_internal[row, col] &&
            Yao.apply!(register,
                Yao.repeat(wire_count, Yao.X, (horizontal_wire,)))
        row > 1 && errors.vertical_internal[row - 1, col] &&
            Yao.apply!(register,
                Yao.repeat(wire_count, Yao.X, (vertical_wire,)))
        col == cols && errors.east_boundary[row] &&
            Yao.apply!(register,
                Yao.repeat(wire_count, Yao.X, (horizontal_wire,)))
        row == 1 && errors.north_boundary[col] &&
            Yao.apply!(register,
                Yao.repeat(wire_count, Yao.X, (vertical_wire,)))
    end

    distribution = Dict{Tuple,Float64}()
    for (basis_state, probability) in zip(Yao.basis(register), Yao.probs(register))
        probability <= 1e-14 && continue
        key = Tuple(Bool(basis_state[wire]) for wire in 1:physical_count)
        distribution[key] = get(distribution, key, 0.0) + probability
    end
    return distribution
end

function distributions_are_approximately_equal(first, second; atol=1e-12)
    keys(first) == keys(second) || return false
    return all(isapprox(first[key], second[key]; atol=atol, rtol=atol)
               for key in keys(first))
end

function with_single_boundary_error(errors, group, index)
    west = falses(length(errors.west_boundary))
    east = falses(length(errors.east_boundary))
    south = falses(length(errors.south_boundary))
    north = falses(length(errors.north_boundary))
    getproperty((; west_boundary=west, east_boundary=east,
                  south_boundary=south, north_boundary=north), group)[index] = true
    return VirtualBondErrors(
        copy(errors.horizontal_internal), copy(errors.vertical_internal),
        west, east, south, north)
end

@testset "Exact small-patch trajectory distributions" begin
    cases = (
        (ToricCodeTrajectoryModel(1, 1), test_virtual_errors(1, 1)),
        (ToricCodeTrajectoryModel(1, 2), test_virtual_errors(
            1, 2; horizontal=BitMatrix(reshape(Bool[1], 1, 1)))),
        (ToricCodeTrajectoryModel(2, 1), test_virtual_errors(
            2, 1; vertical=BitMatrix(reshape(Bool[1], 1, 1)))),
    )
    for (model, errors) in cases
        classical = exact_classical_distribution(model, errors)
        yao = exact_yao_distribution(model, errors)
        @test distributions_are_approximately_equal(classical, yao)
        @test sum(values(classical)) ≈ 1.0
        @test sum(values(yao)) ≈ 1.0
        for group in (
                :west_boundary, :east_boundary,
                :south_boundary, :north_boundary)
            for index in eachindex(getproperty(errors, group))
                boundary_error = with_single_boundary_error(
                    errors, group, index)
                @test distributions_are_approximately_equal(
                    yao, exact_yao_distribution(model, boundary_error))
            end
        end
    end
end

@testset "Virtual-bond trajectory model" begin
    model = ToricCodeTrajectoryModel(3, 4)
    @test model.rows == 3
    @test model.cols == 4
    @test model.layers == [
        [(3, 1)],
        [(3, 2), (2, 1)],
        [(3, 3), (2, 2), (1, 1)],
        [(3, 4), (2, 3), (1, 2)],
        [(2, 4), (1, 3)],
        [(1, 4)],
    ]

    no_errors = sample_virtual_errors(
        MersenneTwister(11), model; error_rate=0.0)
    @test no_errors isa VirtualBondErrors
    @test size(no_errors.horizontal_internal) == (3, 3)
    @test size(no_errors.vertical_internal) == (2, 4)
    @test length(no_errors.west_boundary) == 3
    @test length(no_errors.east_boundary) == 3
    @test length(no_errors.south_boundary) == 4
    @test length(no_errors.north_boundary) == 4
    @test count(no_errors) == 0
    @test length(no_errors) == 31

    all_errors = sample_virtual_errors(
        MersenneTwister(22), model; error_rate=1.0)
    @test count(all_errors) == length(all_errors) == 31
    @test length(all_errors) ==
          model.rows * (model.cols - 1) +
          (model.rows - 1) * model.cols + 2 * model.rows + 2 * model.cols

    probability = 0.31
    draws = 5_000
    selected_bits = falses(draws, 6)
    statistical_rng = MersenneTwister(23)
    for draw in 1:draws
        sampled = sample_virtual_errors(
            statistical_rng, model; error_rate=probability)
        selected_bits[draw, :] .= (
            sampled.horizontal_internal[1, 1],
            sampled.vertical_internal[1, 1],
            sampled.west_boundary[1], sampled.east_boundary[1],
            sampled.south_boundary[1], sampled.north_boundary[1])
    end
    @test all(abs(mean(selected_bits[:, group]) - probability) < 0.025
              for group in 1:6)
    @test all(abs(mean(selected_bits[:, first] .& selected_bits[:, second]) -
                  probability^2) < 0.025
              for first in 1:5 for second in (first + 1):6)

    @test log2_postselection_probability(model) == -7
    @test postselection_probability(model) == 1 / 128

    @test_throws ArgumentError ToricCodeTrajectoryModel(0, 2)
    @test_throws ArgumentError ToricCodeTrajectoryModel(2, 0)
    @test_throws ArgumentError sample_virtual_errors(
        MersenneTwister(1), model; error_rate=-0.1)
    @test_throws ArgumentError sample_virtual_errors(
        MersenneTwister(1), model; error_rate=1.1)
    @test_throws ArgumentError sample_virtual_errors(
        MersenneTwister(1), model; error_rate=NaN)

    valid = test_virtual_errors(3, 4)
    malformed = (
        VirtualBondErrors(
            falses(3, 4), valid.vertical_internal,
            valid.west_boundary, valid.east_boundary,
            valid.south_boundary, valid.north_boundary),
        VirtualBondErrors(
            valid.horizontal_internal, falses(3, 4),
            valid.west_boundary, valid.east_boundary,
            valid.south_boundary, valid.north_boundary),
        VirtualBondErrors(
            valid.horizontal_internal, valid.vertical_internal,
            falses(4), valid.east_boundary,
            valid.south_boundary, valid.north_boundary),
        VirtualBondErrors(
            valid.horizontal_internal, valid.vertical_internal,
            valid.west_boundary, falses(4),
            valid.south_boundary, valid.north_boundary),
        VirtualBondErrors(
            valid.horizontal_internal, valid.vertical_internal,
            valid.west_boundary, valid.east_boundary,
            falses(3), valid.north_boundary),
        VirtualBondErrors(
            valid.horizontal_internal, valid.vertical_internal,
            valid.west_boundary, valid.east_boundary,
            valid.south_boundary, falses(3)),
    )
    for errors in malformed
        @test_throws DimensionMismatch sample_trajectory(
            MersenneTwister(1), model, errors)
    end
end

@testset "Literal Yao measured trajectories" begin
    model = ToricCodeTrajectoryModel(2, 2)
    errors = test_virtual_errors(
        2, 2;
        horizontal=BitMatrix(reshape(Bool[1, 0], 2, 1)),
        vertical=BitMatrix([0 1]),
        west=BitVector([1, 0]), east=BitVector([0, 1]),
        south=BitVector([1, 0]), north=BitVector([0, 1]))
    trajectory = sample_yao_trajectory(
        MersenneTwister(105), model, errors)

    @test trajectory isa ToricCodeTrajectory
    @test trajectory.errors === errors
    @test size(trajectory.measurements) == (2, 2, 4)
    @test all(iseven(sum(trajectory.measurements[row, col, :]))
              for row in 1:2 for col in 1:2)
    @test BitMatrix([
        xor(trajectory.measurements[row, 1, 1],
            trajectory.measurements[row, 2, 3])
        for row in 1:2, _ in 1:1
    ]) == errors.horizontal_internal
    @test BitMatrix([
        xor(trajectory.measurements[2, col, 2],
            trajectory.measurements[1, col, 4])
        for _ in 1:1, col in 1:2
    ]) == errors.vertical_internal

    @test_throws ArgumentError sample_yao_trajectory(
        MersenneTwister(1), model, errors; max_qubits=7)
    @test_throws ArgumentError sample_yao_trajectory(
        MersenneTwister(1), model, errors; max_qubits=0)
end

@testset "Scalable measured trajectories" begin
    model = ToricCodeTrajectoryModel(3, 3)
    horizontal = falses(3, 2)
    horizontal[1, 1] = true
    horizontal[3, 2] = true
    vertical = falses(2, 3)
    vertical[1, 2] = true
    vertical[2, 3] = true
    errors = test_virtual_errors(
        3, 3; horizontal=horizontal, vertical=vertical)

    trajectory = sample_trajectory(MersenneTwister(31), model, errors)
    @test trajectory isa ToricCodeTrajectory
    @test size(trajectory.measurements) == (3, 3, 4)
    @test trajectory.errors === errors
    @test all(iseven(sum(trajectory.measurements[row, col, :]))
              for row in 1:3 for col in 1:3)

    measured_horizontal = BitMatrix([
        xor(trajectory.measurements[row, col, 1],
            trajectory.measurements[row, col + 1, 3])
        for row in 1:3, col in 1:2
    ])
    measured_vertical = BitMatrix([
        xor(trajectory.measurements[row + 1, col, 2],
            trajectory.measurements[row, col, 4])
        for row in 1:2, col in 1:3
    ])
    @test measured_horizontal == horizontal
    @test measured_vertical == vertical

    boundary_errors = test_virtual_errors(
        3, 3;
        horizontal=horizontal,
        vertical=vertical,
        west=trues(3), east=trues(3),
        south=trues(3), north=trues(3))
    same_seed = MersenneTwister(77)
    without_boundary = sample_trajectory(same_seed, model, errors)
    with_boundary = sample_trajectory(
        MersenneTwister(77), model, boundary_errors)
    @test without_boundary.measurements == with_boundary.measurements
    @test without_boundary.errors != with_boundary.errors

    sampled = sample_trajectory(
        MersenneTwister(91), model; error_rate=1.0)
    @test count(sampled.errors) == length(sampled.errors)

    malformed = test_virtual_errors(2, 3)
    @test_throws DimensionMismatch sample_trajectory(
        MersenneTwister(1), model, malformed)
end

@testset "Trajectory mismatch observables" begin
    model = ToricCodeTrajectoryModel(2, 2)
    one_horizontal = test_virtual_errors(
        2, 2;
        horizontal=BitMatrix(reshape(Bool[1, 0], 2, 1)))
    trajectory = sample_trajectory(
        MersenneTwister(151), model, one_horizontal)
    mismatches = bond_mismatches(trajectory)
    @test mismatches.horizontal == one_horizontal.horizontal_internal
    @test mismatches.vertical == one_horizontal.vertical_internal

    observables = trajectory_observables(trajectory)
    @test observables isa TrajectoryObservables
    @test observables.sampled_error_density == 1 / 12
    @test observables.boundary_error_density == 0.0
    @test observables.mismatch_density == 1 / 4
    @test observables.frustration_density == 1.0
    @test observables.largest_cluster_fraction == 1 / 2
    @test observables.spans_horizontal
    @test !observables.spans_vertical

    boundary_only = test_virtual_errors(
        2, 2;
        west=trues(2), east=trues(2),
        south=trues(2), north=trues(2))
    boundary_observables = trajectory_observables(sample_trajectory(
        MersenneTwister(152), model, boundary_only))
    @test boundary_observables.sampled_error_density == 2 / 3
    @test boundary_observables.boundary_error_density == 1.0
    @test boundary_observables.mismatch_density == 0.0
    @test boundary_observables.frustration_density == 0.0
    @test boundary_observables.largest_cluster_fraction == 1 / 4
    @test !boundary_observables.spans_horizontal
    @test !boundary_observables.spans_vertical

    all_errors = sample_virtual_errors(
        MersenneTwister(153), model; error_rate=1.0)
    all_observables = trajectory_observables(sample_trajectory(
        MersenneTwister(154), model, all_errors))
    @test all_observables.sampled_error_density == 1.0
    @test all_observables.boundary_error_density == 1.0
    @test all_observables.mismatch_density == 1.0
    @test all_observables.frustration_density == 0.0
    @test all_observables.largest_cluster_fraction == 1.0
    @test all_observables.spans_horizontal
    @test all_observables.spans_vertical

    singleton = ToricCodeTrajectoryModel(1, 1)
    singleton_observables = trajectory_observables(sample_trajectory(
        MersenneTwister(155), singleton; error_rate=0.0))
    @test singleton_observables.mismatch_density == 0.0
    @test singleton_observables.frustration_density == 0.0
    @test singleton_observables.largest_cluster_fraction == 1.0
    @test !singleton_observables.spans_horizontal
    @test !singleton_observables.spans_vertical
end
