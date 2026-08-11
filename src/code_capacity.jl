"""Open rectangular patch used for data-edge code-capacity experiments."""
struct OpenCodeCapacityModel
    rows::Int
    cols::Int

    function OpenCodeCapacityModel(rows::Integer, cols::Integer)
        rows >= 2 || throw(ArgumentError("rows must be at least 2"))
        cols >= 2 || throw(ArgumentError("cols must be at least 2"))
        new(Int(rows), Int(cols))
    end
end

"""Independent data-edge errors on an open rectangular patch."""
struct DataEdgeErrors
    horizontal::BitMatrix
    vertical::BitMatrix

    function DataEdgeErrors(horizontal::BitMatrix, vertical::BitMatrix)
        rows, horizontal_cols = size(horizontal)
        vertical_rows, cols = size(vertical)
        rows >= 2 || throw(ArgumentError("data-edge arrays must have at least 2 rows"))
        cols >= 2 || throw(ArgumentError("data-edge arrays must have at least 2 columns"))
        horizontal_cols == cols - 1 ||
            throw(ArgumentError("horizontal errors must have size (R, C - 1)"))
        vertical_rows == rows - 1 ||
            throw(ArgumentError("vertical errors must have size (R - 1, C)"))
        new(horizontal, vertical)
    end
end

"""Sample independent data-edge errors at a common Bernoulli error rate."""
function sample_data_edge_errors(rng, model::OpenCodeCapacityModel; error_rate::Real)
    0 <= error_rate <= 1 || throw(ArgumentError("error_rate must lie in [0, 1]"))
    horizontal = BitMatrix(rand(rng, model.rows, model.cols - 1) .< error_rate)
    vertical = BitMatrix(rand(rng, model.rows - 1, model.cols) .< error_rate)
    return DataEdgeErrors(horizontal, vertical)
end

"""Return the plaquette syndrome induced by independent data-edge errors."""
function code_capacity_syndrome(errors::DataEdgeErrors)
    rows, cols = size(errors.horizontal, 1), size(errors.vertical, 2)
    syndrome = falses(rows - 1, cols - 1)
    for row in 1:(rows - 1), col in 1:(cols - 1)
        syndrome[row, col] = xor(errors.horizontal[row, col], errors.horizontal[row + 1, col],
                                 errors.vertical[row, col], errors.vertical[row, col + 1])
    end
    return syndrome
end

"""Return a central transverse data-edge mask for the requested logical sector."""
function logical_cut(model::OpenCodeCapacityModel; sector::Symbol)
    horizontal = falses(model.rows, model.cols - 1)
    vertical = falses(model.rows - 1, model.cols)
    if sector === :north_south
        vertical[cld(model.rows - 1, 2), :] .= true
    elseif sector === :east_west
        horizontal[:, cld(model.cols - 1, 2)] .= true
    else
        throw(ArgumentError("sector must be :north_south or :east_west"))
    end
    return (horizontal=horizontal, vertical=vertical)
end
