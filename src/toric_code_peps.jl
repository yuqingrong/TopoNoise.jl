const _DIRECTIONS = (:east, :north, :west, :south)

"""
    ToricCodePEPS

A finite rectangular PEPS for the doubled-edge toric-code state. `tensors` is
the matrix of local `ITensor`s and `physical_indices[row, col]` stores the four
physical qubit indices in `(east, north, west, south)` order.
"""
struct ToricCodePEPS
    tensors::Matrix{ITensor}
    physical_indices::Matrix{NTuple{4,Index}}
end

Base.size(peps::ToricCodePEPS) = size(peps.tensors)
Base.size(peps::ToricCodePEPS, dim::Integer) = size(peps.tensors, dim)
Base.axes(peps::ToricCodePEPS) = axes(peps.tensors)
Base.getindex(peps::ToricCodePEPS, row::Integer, col::Integer) = peps.tensors[row, col]

"""
    physicalinds(peps, row, col)

Return the four physical qubit indices at `(row, col)` in
`(east, north, west, south)` order.
"""
function physicalinds(peps::ToricCodePEPS, row::Integer, col::Integer)
    checkbounds(peps.physical_indices, row, col)
    return peps.physical_indices[row, col]
end

raw"""
    toric_code_local_tensor([T=Float64])

Return the rank-eight local tensor

```math
T^{i j k l}_{\alpha \beta \gamma \delta}
= \frac{1}{\sqrt{2}}
  \delta_{i\alpha}\delta_{j\beta}\delta_{k\gamma}\delta_{l\delta}
```

when `\alpha + \beta + \gamma + \delta` is even, and zero otherwise. The
index order is the four physical indices `(E, N, W, S)`, followed by the four
virtual indices `(α, β, γ, δ)`. Binary values 0 and 1 occupy Julia index
positions 1 and 2, respectively.
"""
function toric_code_local_tensor(::Type{T}=Float64) where {T<:Number}
    _supports_normalization(T) || throw(ArgumentError(
        "T must be a real or complex floating-point type, got $T"))
    tensor = zeros(T, ntuple(_ -> 2, 8))
    coefficient = inv(sqrt(convert(T, 2)))
    for east in 0:1, north in 0:1, west in 0:1, south in 0:1
        iseven(east + north + west + south) || continue
        bits = (east + 1, north + 1, west + 1, south + 1)
        tensor[bits..., bits...] = coefficient
    end
    return tensor
end

function _physical_indices(rows::Int, cols::Int)
    physical = Matrix{NTuple{4,Index}}(undef, rows, cols)
    for row in 1:rows, col in 1:cols
        physical[row, col] = ntuple(4) do direction
            name = _DIRECTIONS[direction]
            Index(2, "Site,$name,r=$row,c=$col")
        end
    end
    return physical
end

function _virtual_indices(rows::Int, cols::Int)
    vertical = Matrix{Index}(undef, rows - 1, cols)
    for row in axes(vertical, 1), col in axes(vertical, 2)
        vertical[row, col] = Index(2, "Link,v,r=$row,c=$col")
    end

    horizontal = Matrix{Index}(undef, rows, cols - 1)
    for row in axes(horizontal, 1), col in axes(horizontal, 2)
        horizontal[row, col] = Index(2, "Link,h,r=$row,c=$col")
    end
    return vertical, horizontal
end

function _boundary_index(row::Int, col::Int, direction::Symbol)
    return Index(2, "Boundary,$direction,r=$row,c=$col")
end

function _site_virtual_indices(
        row::Int, col::Int, rows::Int, cols::Int,
        vertical::Matrix{Index}, horizontal::Matrix{Index})
    north = row > 1 ? vertical[row - 1, col] : _boundary_index(row, col, :north)
    east = col < cols ? horizontal[row, col] : _boundary_index(row, col, :east)
    south = row < rows ? vertical[row, col] : _boundary_index(row, col, :south)
    west = col > 1 ? horizontal[row, col - 1] : _boundary_index(row, col, :west)
    virtual = (east, north, west, south)
    boundary = (col == cols, row == 1, col == 1, row == rows)
    return virtual, boundary
end

function _x_boundary_vector(::Type{T}, index::Index) where {T<:Number}
    plus = ITensor(T, index)
    coefficient = inv(sqrt(convert(T, 2)))
    plus[index => 1] = coefficient
    plus[index => 2] = coefficient
    return plus
end

function _supports_normalization(::Type{T}) where {T}
    return T <: AbstractFloat || T <: Complex{<:AbstractFloat}
end

"""
    toric_code_peps(rows, cols; scalar_type=Float64)

Construct a normalized finite `rows × cols` doubled-edge toric-code PEPS.
Each site has four physical qubits ordered `(east, north, west, south)`.
Internal virtual indices are shared by neighboring tensors, which copies an
internal bond bit into the two physical states `|i i⟩`. Each dangling virtual
index is contracted with the x-polarized state `(|0⟩ + |1⟩)/√2`.

The returned network is normalized analytically. If `V = rows*cols`, `E` is
the number of internal edges, and `B = 2rows + 2cols`, its
`2^(E + B - V)` allowed basis configurations all have the same amplitude.
"""
function toric_code_peps(
        rows::Integer, cols::Integer; scalar_type::Type{T}=Float64) where {T<:Number}
    rows > 0 || throw(ArgumentError("rows must be positive, got $rows"))
    cols > 0 || throw(ArgumentError("cols must be positive, got $cols"))
    _supports_normalization(T) || throw(ArgumentError(
        "scalar_type must be a real or complex floating-point type, got $T"))

    nrows, ncols = Int(rows), Int(cols)
    physical = _physical_indices(nrows, ncols)
    vertical, horizontal = _virtual_indices(nrows, ncols)
    local_data = toric_code_local_tensor(T)
    tensors = Matrix{ITensor}(undef, nrows, ncols)

    vertices = nrows * ncols
    internal_edges = nrows * (ncols - 1) + (nrows - 1) * ncols
    site_scale = convert(T, exp2(1 - internal_edges / (2 * vertices)))

    for row in 1:nrows, col in 1:ncols
        virtual, boundary = _site_virtual_indices(
            row, col, nrows, ncols, vertical, horizontal)
        tensor = ITensor(T, local_data, physical[row, col]..., virtual...)
        for direction in 1:4
            boundary[direction] || continue
            tensor *= _x_boundary_vector(T, virtual[direction])
        end
        tensors[row, col] = site_scale * tensor
    end

    return ToricCodePEPS(tensors, physical)
end
