"""
Experimental sparse zip-up microkernel for one compiled Fourier row branch.

This file deliberately has no decoder integration.  It converts an existing
`FourierRowMPO` into compact local transition tables and applies one fixed
seam branch to one dense-array MPS.  It exists only to measure whether a
custom local kernel can pass the row-level performance gate before replacing
any of the reference or production code.
"""
module TNMLRowStripSparseMicrokernel

using LinearAlgebra, ITensors, ITensorMPS
using ..TNMLRowStripMPS

export DenseMPSChain, SparseRowOperator, SparseZipupWorkspace,
       SparseZipupTimings,
       compile_sparse_row_operator, dense_mps_chain,
       sparse_zipup_apply!, dense_chain_vector

struct DenseMPSChain
    tensors::Vector{Array{Float64,3}}
end

struct SparseSiteOperator
    left_dimension::Int
    input_dimension::Int
    output_dimension::Int
    right_dimension::Int
    offsets::Vector{Int}
    right_states::Vector{Int}
    output_states::Vector{Int}
    weights::Vector{Float64}
end

struct SparseRowOperator
    sites::Vector{SparseSiteOperator}
    maximum_transitions_per_input::Int
    transition_count::Int
end

mutable struct SparseZipupWorkspace
    matrix::Matrix{Float64}
    carry::Array{Float64,3}
    output_tensors::Vector{Array{Float64,3}}
    bond_upper_bounds::Vector{Int}
    actual_bonds::Vector{Int}
    maxdim::Int
end

struct SparseZipupResult
    workspace::SparseZipupWorkspace
end

mutable struct SparseZipupTimings
    accumulation_ns::Int
    factorization_ns::Int
    compression_ns::Int
    final_contraction_ns::Int
end

SparseZipupTimings() = SparseZipupTimings(0, 0, 0, 0)

function _reset!(timings::SparseZipupTimings)
    timings.accumulation_ns = 0
    timings.factorization_ns = 0
    timings.compression_ns = 0
    timings.final_contraction_ns = 0
    return timings
end

function _ordered_tensor_array(tensor::ITensor, left, physical, right)
    ordered = Index[]
    isnothing(left) || push!(ordered, left)
    push!(ordered, physical)
    isnothing(right) || push!(ordered, right)
    values = Array(tensor, ordered...)
    return reshape(Float64.(values),
        isnothing(left) ? 1 : dim(left), dim(physical),
        isnothing(right) ? 1 : dim(right))
end

"""Copy an ITensor MPS into `(left, physical, right)` dense arrays."""
function dense_mps_chain(state::ITensorMPS.MPS)
    count = length(state)
    tensors = Array{Float64,3}[]
    sizehint!(tensors, count)
    for position in 1:count
        left = position == 1 ? nothing : ITensorMPS.linkind(state, position - 1)
        right = position == count ? nothing : ITensorMPS.linkind(state, position)
        physical = only(ITensorMPS.siteinds(state, position))
        push!(tensors,
            _ordered_tensor_array(state[position], left, physical, right))
    end
    return DenseMPSChain(tensors)
end

function _operator_array(row, state, position)
    count = length(row.operator)
    tensor = row.operator[position]
    left = position == 1 ? nothing :
        ITensorMPS.linkind(row.operator, position - 1)
    right = position == count ? nothing :
        ITensorMPS.linkind(row.operator, position)
    common = collect(commoninds(tensor, state[position]))
    length(common) == 1 || error(
        "expected exactly one MPO/MPS physical input index")
    input = only(common)
    output = prime(input)
    hasind(tensor, output) || error("MPO output index is not prime(input)")
    ordered = Index[]
    isnothing(left) || push!(ordered, left)
    push!(ordered, input)
    push!(ordered, output)
    isnothing(right) || push!(ordered, right)
    values = Array(tensor, ordered...)
    return reshape(Float64.(values),
        isnothing(left) ? 1 : dim(left), dim(input), dim(output),
        isnothing(right) ? 1 : dim(right))
end

function _compile_site(values::Array{Float64,4})
    left_dimension, input_dimension, output_dimension, right_dimension =
        size(values)
    offsets = Int[1]
    right_states = Int[]
    output_states = Int[]
    weights = Float64[]
    maximum = 0
    for input_state in 1:input_dimension, left_state in 1:left_dimension
        before = length(weights)
        for right_state in 1:right_dimension,
                output_state in 1:output_dimension
            weight = values[left_state, input_state, output_state, right_state]
            iszero(weight) && continue
            push!(right_states, right_state)
            push!(output_states, output_state)
            push!(weights, weight)
        end
        push!(offsets, length(weights) + 1)
        maximum = max(maximum, length(weights) - before)
    end
    site = SparseSiteOperator(left_dimension, input_dimension,
        output_dimension, right_dimension, offsets, right_states,
        output_states, weights)
    return site, maximum
end

"""Compile one raw finite-state MPO into CSR-like local transitions."""
function compile_sparse_row_operator(row::TNMLRowStripMPS.FourierRowMPO,
        state::ITensorMPS.MPS)
    length(row.operator) == length(state) || throw(DimensionMismatch(
        "row MPO and input MPS have different lengths"))
    sites = SparseSiteOperator[]
    maximum = 0
    transitions = 0
    for position in 1:length(state)
        site, local_maximum = _compile_site(
            _operator_array(row, state, position))
        push!(sites, site)
        maximum = max(maximum, local_maximum)
        transitions += length(site.weights)
    end
    return SparseRowOperator(sites, maximum, transitions)
end

function _validate_chain(operator::SparseRowOperator, input::DenseMPSChain)
    length(operator.sites) == length(input.tensors) ||
        throw(DimensionMismatch("operator and MPS have different lengths"))
    for position in eachindex(operator.sites)
        site = operator.sites[position]
        tensor = input.tensors[position]
        size(tensor, 2) == site.input_dimension || throw(DimensionMismatch(
            "physical input dimension differs at site $position"))
        position == 1 || size(tensor, 1) ==
            size(input.tensors[position - 1], 3) || throw(DimensionMismatch(
                "MPS bond dimension differs at site $position"))
        position == 1 || site.left_dimension ==
            operator.sites[position - 1].right_dimension ||
            throw(DimensionMismatch(
                "MPO bond dimension differs at site $position"))
    end
    size(first(input.tensors), 1) == 1 || error("left MPS boundary is open")
    first(operator.sites).left_dimension == 1 ||
        error("left MPO boundary is open")
    size(last(input.tensors), 3) == 1 || error("right MPS boundary is open")
    last(operator.sites).right_dimension == 1 ||
        error("right MPO boundary is open")
    return nothing
end

function SparseZipupWorkspace(operator::SparseRowOperator,
        input::DenseMPSChain; maxdim::Integer)
    maxdim > 0 || throw(ArgumentError("maxdim must be positive"))
    _validate_chain(operator, input)
    count = length(operator.sites)
    count >= 2 || throw(ArgumentError(
        "the row microkernel expects at least two sites"))
    bounds = ones(Int, count - 1)
    left_rank = 1
    maximum_rows = 1
    maximum_columns = 1
    maximum_input_left = 1
    maximum_operator_left = 1
    for position in 1:count-1
        site = operator.sites[position]
        tensor = input.tensors[position]
        rows = Base.checked_mul(left_rank, site.output_dimension)
        columns = Base.checked_mul(size(tensor, 3), site.right_dimension)
        right_rank = min(Int(maxdim), rows, columns)
        bounds[position] = right_rank
        maximum_rows = max(maximum_rows, rows)
        maximum_columns = max(maximum_columns, columns)
        maximum_input_left = max(maximum_input_left, size(tensor, 1))
        maximum_operator_left = max(maximum_operator_left,
            site.left_dimension)
        left_rank = right_rank
    end
    final_site = last(operator.sites)
    final_tensor = last(input.tensors)
    maximum_rows = max(maximum_rows,
        Base.checked_mul(left_rank, final_site.output_dimension))
    maximum_input_left = max(maximum_input_left, size(final_tensor, 1))
    maximum_operator_left = max(maximum_operator_left,
        final_site.left_dimension)

    output_tensors = Array{Float64,3}[]
    for position in 1:count
        left = position == 1 ? 1 : bounds[position - 1]
        right = position == count ? 1 : bounds[position]
        push!(output_tensors, zeros(Float64, left,
            operator.sites[position].output_dimension, right))
    end
    carry_rank = max(1, maximum(bounds; init=1))
    return SparseZipupWorkspace(
        zeros(Float64, maximum_rows, maximum_columns),
        zeros(Float64, carry_rank, maximum_input_left,
            maximum_operator_left),
        output_tensors, bounds, ones(Int, count - 1), Int(maxdim))
end

@inline function _transition_range(site::SparseSiteOperator,
        left_state::Int, input_state::Int)
    key = left_state + (input_state - 1) * site.left_dimension
    return site.offsets[key]:(site.offsets[key + 1] - 1)
end

function _accumulate_site!(matrix, carry, left_rank, input_tensor,
        site::SparseSiteOperator)
    input_left, input_dimension, input_right = size(input_tensor)
    rows = left_rank * site.output_dimension
    columns = input_right * site.right_dimension
    target = @view matrix[1:rows, 1:columns]
    fill!(target, 0.0)
    for input_state in 1:input_dimension,
            left_state in 1:site.left_dimension
        left_block = @view carry[1:left_rank, 1:input_left, left_state]
        input_block = @view input_tensor[:, input_state, :]
        for transition in _transition_range(site, left_state, input_state)
            right_state = site.right_states[transition]
            output_state = site.output_states[transition]
            row_range = (output_state - 1) * left_rank .+ (1:left_rank)
            column_range = (right_state - 1) * input_right .+ (1:input_right)
            destination = @view target[row_range, column_range]
            mul!(destination, left_block, input_block,
                site.weights[transition], 1.0)
        end
    end
    return target
end

function _retained_rank(singular_values, maxdim::Int, cutoff::Float64)
    rank = min(maxdim, length(singular_values))
    cutoff <= 0 && return rank
    total = sum(abs2, singular_values)
    iszero(total) && return 1
    discarded = sum(abs2, @view singular_values[rank + 1:end])
    while rank > 1
        candidate = discarded + abs2(singular_values[rank])
        candidate <= cutoff * total || break
        discarded = candidate
        rank -= 1
    end
    return rank
end

"""
Apply one sparse fixed-seam branch with a left-to-right zip-up sweep.

The workspace owns every large accumulation and output buffer.  LAPACK still
owns the temporary factorization buffers returned by `svd!`; the benchmark
records those allocations explicitly.
"""
function sparse_zipup_apply!(workspace::SparseZipupWorkspace,
        operator::SparseRowOperator, input::DenseMPSChain;
        maxdim::Integer=workspace.maxdim, cutoff::Real=0.0,
        timings::Union{Nothing,SparseZipupTimings}=nothing)
    0 < maxdim <= workspace.maxdim || throw(ArgumentError(
        "maxdim exceeds the workspace capacity"))
    isfinite(cutoff) && 0 <= cutoff < 1 || throw(ArgumentError(
        "cutoff must be finite and in [0,1)"))
    _validate_chain(operator, input)
    isnothing(timings) || _reset!(timings)
    fill!(workspace.carry, 0.0)
    workspace.carry[1, 1, 1] = 1.0
    left_rank = 1
    count = length(operator.sites)

    for position in 1:count-1
        site = operator.sites[position]
        start = time_ns()
        matrix = _accumulate_site!(workspace.matrix, workspace.carry,
            left_rank, input.tensors[position], site)
        isnothing(timings) ||
            (timings.accumulation_ns += time_ns() - start)
        start = time_ns()
        factorization = svd!(matrix; full=false)
        right_rank = _retained_rank(
            factorization.S, Int(maxdim), Float64(cutoff))
        isnothing(timings) ||
            (timings.factorization_ns += time_ns() - start)
        start = time_ns()
        workspace.actual_bonds[position] = right_rank
        output = workspace.output_tensors[position]
        fill!(output, 0.0)
        @views reshape(output[1:left_rank, :, 1:right_rank],
            left_rank * site.output_dimension, right_rank) .=
            factorization.U[:, 1:right_rank]

        input_right = size(input.tensors[position], 3)
        fill!(workspace.carry, 0.0)
        @inbounds for right_state in 1:site.right_dimension,
                input_state in 1:input_right,
                rank in 1:right_rank
            column = input_state + (right_state - 1) * input_right
            workspace.carry[rank, input_state, right_state] =
                factorization.S[rank] * factorization.Vt[rank, column]
        end
        isnothing(timings) ||
            (timings.compression_ns += time_ns() - start)
        left_rank = right_rank
    end

    site = last(operator.sites)
    start = time_ns()
    matrix = _accumulate_site!(workspace.matrix, workspace.carry,
        left_rank, last(input.tensors), site)
    size(matrix, 2) == 1 || error("right boundary did not close")
    output = last(workspace.output_tensors)
    fill!(output, 0.0)
    @views reshape(output[1:left_rank, :, 1], :) .= matrix[:, 1]
    isnothing(timings) ||
        (timings.final_contraction_ns += time_ns() - start)
    return SparseZipupResult(workspace)
end

function _result_tensor(result::SparseZipupResult, position::Int)
    workspace = result.workspace
    count = length(workspace.output_tensors)
    left = position == 1 ? 1 : workspace.actual_bonds[position - 1]
    right = position == count ? 1 : workspace.actual_bonds[position]
    return @view workspace.output_tensors[position][1:left, :, 1:right]
end

"""Materialize the physical coefficient vector of a sparse zip-up result."""
function dense_chain_vector(result::SparseZipupResult)
    count = length(result.workspace.output_tensors)
    first_tensor = _result_tensor(result, 1)
    coefficients = reshape(copy(first_tensor),
        size(first_tensor, 2), size(first_tensor, 3))
    physical_dimension = size(coefficients, 1)
    for position in 2:count
        tensor = _result_tensor(result, position)
        size(coefficients, 2) == size(tensor, 1) || error(
            "result MPS bond mismatch at site $position")
        product = coefficients * reshape(tensor,
            size(tensor, 1), size(tensor, 2) * size(tensor, 3))
        physical_dimension *= size(tensor, 2)
        coefficients = reshape(product,
            physical_dimension, size(tensor, 3))
    end
    size(coefficients, 2) == 1 || error("result MPS did not close")
    return vec(coefficients)
end

end

