function _validate_peps_graph_input(peps::ToricCodePEPS)
    size(peps.tensors) == size(peps.physical_indices) || throw(DimensionMismatch(
        "tensor grid has size $(size(peps.tensors)), but physical-index grid " *
        "has size $(size(peps.physical_indices))"))

    rows, cols = size(peps)
    rows > 0 && cols > 0 || throw(ArgumentError(
        "the PEPS topology must be a nonempty rectangular grid"))
    index_sites = Dict{Index,Vector{Tuple{Int,Int}}}()
    physical_sets = Dict{Tuple{Int,Int},Set{Index}}()

    for row in 1:rows, col in 1:cols
        site = (row, col)
        tensor = peps[row, col]
        physical = physicalinds(peps, row, col)
        length(unique(physical)) == 4 || throw(ArgumentError(
            "site $site must have four distinct physical indices"))

        for (direction, index) in zip(_DIRECTIONS, physical)
            ITensors.hastags(index, "Site,$direction") || throw(ArgumentError(
                "physical index $index at site $site is not tagged for $direction"))
            ITensors.hastags(index, "r=$row,c=$col") || throw(ArgumentError(
                "physical index $index has coordinates inconsistent with site $site"))
            index in ITensors.inds(tensor) || throw(ArgumentError(
                "physical index $index at site $site is absent from its tensor"))
        end

        physical_sets[site] = Set(physical)
        for index in ITensors.inds(tensor)
            ITensors.dim(index) == 2 || throw(ArgumentError(
                "index $index at site $site has dimension $(ITensors.dim(index)); " *
                "toric-code PEPS graphs require qubit indices of dimension 2"))
            push!(get!(index_sites, index, Tuple{Int,Int}[]), site)
        end
    end

    for (index, sites) in index_sites
        length(sites) <= 2 || throw(ArgumentError(
            "index $index is shared by $(length(sites)) tensors; hyperedges are unsupported"))
        if length(sites) == 1
            site = only(sites)
            index in physical_sets[site] || throw(ArgumentError(
                "unexpected external non-physical index $index at site $site"))
        else
            first_site, second_site = sites
            distance = abs(first_site[1] - second_site[1]) +
                       abs(first_site[2] - second_site[2])
            distance == 1 || throw(ArgumentError(
                "index $index connects non-neighboring sites $first_site and $second_site"))
            (index in physical_sets[first_site] || index in physical_sets[second_site]) &&
                throw(ArgumentError(
                    "physical index $index is unexpectedly shared between sites"))
        end
    end

    for row in 1:rows, col in 1:(cols - 1)
        shared = collect(ITensors.commoninds(peps[row, col], peps[row, col + 1]))
        length(shared) == 1 || throw(ArgumentError(
            "horizontal neighbors $((row, col)) and $((row, col + 1)) must share " *
            "exactly one virtual index, found $(length(shared))"))
    end
    for row in 1:(rows - 1), col in 1:cols
        shared = collect(ITensors.commoninds(peps[row, col], peps[row + 1, col]))
        length(shared) == 1 || throw(ArgumentError(
            "vertical neighbors $((row, col)) and $((row + 1, col)) must share " *
            "exactly one virtual index, found $(length(shared))"))
    end
    return nothing
end

"""
    peps_graph(peps::ToricCodePEPS)

Convert a finite toric-code PEPS to an `ITensorNetwork` whose vertices are
labelled by `(row, col)`. Shared virtual indices define nearest-neighbor edges,
while the four directional physical indices remain external legs.
"""
function peps_graph(peps::ToricCodePEPS)
    _validate_peps_graph_input(peps)
    rows, cols = size(peps)
    tensors = Dict{Tuple{Int,Int},ITensor}(
        (row, col) => peps[row, col] for row in 1:rows for col in 1:cols)
    return ITensorNetwork(tensors)
end
