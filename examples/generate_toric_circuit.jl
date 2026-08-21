module GenerateToricCircuitExample

using TopoNoise
import Yao

const DEFAULT_OUTPUT_DIR = normpath(joinpath(
    @__DIR__, "..", "results", "toric-circuit"))

const USAGE =
    "Usage: julia --project=. examples/generate_toric_circuit.jl " *
    "[rows cols] (rows and cols must be positive integers)"

function _dimensions(args)
    isempty(args) && return (2, 2)
    length(args) == 2 ||
        throw(ArgumentError("expected zero or two dimensions"))
    dimensions = try
        parse.(Int, args)
    catch
        throw(ArgumentError("rows and cols must be positive integers"))
    end
    all(>(0), dimensions) ||
        throw(ArgumentError("rows and cols must be positive integers"))
    return Tuple(dimensions)
end

function _site_text(site)
    row, col = site
    return "($row,$col)"
end

function generate(
        rows::Integer, cols::Integer;
        output_dir::AbstractString=DEFAULT_OUTPUT_DIR,
        io::IO=stdout)
    peps = toric_code_peps(rows, cols)
    circuit = sequential_circuit_graph(peps)
    block = yao_circuit(circuit)

    println(io, "lattice: $(rows)x$(cols)")
    println(io, "qubits: $(Yao.nqubits(block))")
    for (number, sites) in enumerate(circuit_layers(circuit))
        println(io, "layer $number: ", join(_site_text.(sites), ", "))
    end

    mkpath(output_dir)
    basename = "toric_circuit_$(rows)x$(cols)"
    outputs = String[]
    for format in (:svg, :pdf, :png)
        path = abspath(joinpath(output_dir, "$basename.$format"))
        Yao.vizcircuit(block; format=format, filename=path)
        push!(outputs, path)
        println(io, "wrote: $path")
    end
    return (; peps, circuit, block, outputs)
end

function main(
        args=ARGS;
        output_dir::AbstractString=DEFAULT_OUTPUT_DIR,
        io::IO=stdout,
        error_io::IO=stderr)
    dimensions = try
        _dimensions(args)
    catch error
        error isa ArgumentError || rethrow()
        println(error_io, "error: ", error.msg)
        println(error_io, USAGE)
        return 1
    end
    generate(dimensions...; output_dir=output_dir, io=io)
    return 0
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(GenerateToricCircuitExample.main())
end
