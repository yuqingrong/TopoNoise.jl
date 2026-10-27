using TopoNoise

function main(args=ARGS)
    if any(arg -> arg in ("-h", "--help"), args)
        println("Usage: julia --project=. examples/toric3d_z_tnml.jl [L=3] [p=0.01]")
        println("Decode one deterministic single-fault syndrome; no sampling or scan.")
        return nothing
    end
    length(args) <= 2 || throw(ArgumentError("expected at most L and p"))
    L = isempty(args) ? 3 : parse(Int, args[1])
    p = length(args) < 2 ? 0.01 : parse(Float64, args[2])

    code = Toric3DEncodedSlice(L, L, L)
    model = build_toric3d_z_fault_model(code)
    topology = build_toric3d_z_tn_topology(model)
    settings = Toric3DZTNSettings()
    decoder = build_toric3d_z_tn_decoder(topology, p; settings)

    # The first raw circuit-fault column gives a reproducible, nonzero syndrome.
    # Decoding conditions on this syndrome and sums over all compatible faults.
    raw_column = 1
    syndrome = BitVector(topology.raw_detector_matrix[:, raw_column] .!= 0)
    println("L = $L, p = $p, maximum bond dimension = $(settings.max_bond_dimension)")
    println("Syndrome: raw circuit-fault column $raw_column; detectors = $(findall(syndrome))")
    flush(stdout)

    result = decode_toric3d_z_tn(decoder, syndrome)
    println("Numerically admissible: ", result.valid)
    println("Approximate posterior in order (0,0), (1,0), (0,1), (1,1):")
    println(result.probabilities)
    if result.valid
        logical_bits = (result.logical_mask & 1, (result.logical_mask >> 1) & 1)
        println("Chosen logical mask: $(result.logical_mask); logical bits: $logical_bits")
    else
        println("No admissible logical choice (logical_mask = -1).")
    end
    println("Approximate log syndrome probability: ", result.log_syndrome)
    println("Maximum bond dimension used: ", result.diagnostics.max_bond_dimension)
    println("The decision is a logical class; valid does not certify convergence.")
    return result
end

main()
