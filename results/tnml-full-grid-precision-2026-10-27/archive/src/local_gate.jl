function _unitary_completion_rtol(::Type{T}) where {T}
    return sqrt(eps(real(float(oneunit(T)))))
end

function _supports_unitary_completion(::Type{T}) where {T}
    return T <: AbstractFloat || T <: Complex{<:AbstractFloat}
end

"""
    unitary_completion(isometry; atol=0, rtol=nothing)

Complete the columns of a tall isometric matrix to a square unitary matrix.
The supplied columns are preserved exactly as the leading columns of the
result. The remaining orthonormal basis is not unique and is chosen from the
nullspace of `adjoint(isometry)`.

The input must use a real or complex floating-point element type, have at
least as many rows as columns, and satisfy `isometry' * isometry ≈ I` under
the requested tolerances. When `rtol` is omitted, it defaults to the square
root of machine epsilon for the input scalar type.
"""
function unitary_completion(
        isometry::AbstractMatrix{T}; atol::Real=0,
        rtol::Union{Nothing,Real}=nothing) where {T}
    _supports_unitary_completion(T) || throw(ArgumentError(
        "isometry must have a real or complex floating-point element type, got $T"))

    rows, columns = size(isometry)
    rows >= columns || throw(DimensionMismatch(
        "isometry must have at least as many rows as columns, got size $(size(isometry))"))

    atol >= 0 || throw(ArgumentError("atol must be nonnegative, got $atol"))
    effective_rtol = isnothing(rtol) ? _unitary_completion_rtol(T) : rtol
    effective_rtol >= 0 || throw(ArgumentError(
        "rtol must be nonnegative, got $effective_rtol"))

    identity_columns = Matrix{T}(I, columns, columns)
    gram = adjoint(isometry) * isometry
    isapprox(gram, identity_columns; atol=atol, rtol=effective_rtol) ||
        throw(ArgumentError(
            "isometry columns must be orthonormal under the requested tolerances"))

    complement = try
        nullspace(adjoint(isometry); atol=atol, rtol=effective_rtol)
    catch error
        error isa MethodError || rethrow()
        throw(ArgumentError(
            "element type $T is not supported by LinearAlgebra.nullspace"))
    end
    expected_columns = rows - columns
    expected_size = (rows, expected_columns)
    size(complement) == expected_size || throw(ArgumentError(
        "expected an orthogonal complement of size $expected_size, " *
        "got $(size(complement)); use tighter tolerances"))

    return hcat(isometry, complement)
end

"""
    toric_code_local_gate([T=Float64]; atol=0, rtol=nothing)

Return a rank-12 unitary completion of [`toric_code_local_tensor`](@ref).
Its indices are ordered as the six outputs
`(i, j, k, l, α, β)`, followed by the six inputs
`(phyᵢ, phyⱼ, phyₖ, phyₗ, γ, δ)`. Every index has dimension two, so the
result represents a six-qubit gate.

Julia index position `1` represents qubit value zero. Consequently, fixing
the four physical input indices to `1` recovers the local tensor exactly:

```julia
gate[:, :, :, :, :, :, 1, 1, 1, 1, :, :] == toric_code_local_tensor(T)
```

The action of the gate outside this physical-zero input subspace depends on
the non-unique orthogonal basis chosen by [`unitary_completion`](@ref).
"""
function toric_code_local_gate(
        ::Type{T}=Float64; atol::Real=0,
        rtol::Union{Nothing,Real}=nothing) where {T<:Number}
    tensor = toric_code_local_tensor(T)
    isometry = reshape(tensor, 64, 4)
    unitary = unitary_completion(isometry; atol=atol, rtol=rtol)

    # The leading four columns vary (γ, δ) while the physical inputs are
    # fixed to zero. Reshape in that order, then expose the diagram order
    # with the four physical inputs before (γ, δ).
    gamma_delta_first = reshape(unitary, ntuple(_ -> 2, 12))
    return permutedims(
        gamma_delta_first, (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 7, 8))
end
