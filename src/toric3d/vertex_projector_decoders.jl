"""Uniform-prior BP+OSD data-X recovery for one declared prepared support."""
struct Toric3DVertexBPOSDDecoder
    window::Toric3DVertexFeedbackWindow
    p::Float64
    p_eff::Float64
    settings::Toric3DDecoderSettings
    bposd::PythonCall.Py
    numpy::PythonCall.Py
    version::String
end
struct Toric3DVertexDecoderFamily
    code::Toric3DVertexProjectorCode
    p::Float64
    settings::Toric3DDecoderSettings
    x_one_step::Vector{Toric3DVertexMatchingDecoder}
    z_one_step::Vector{Toric3DVertexMatchingDecoder}
    x_window::Vector{Toric3DVertexBPOSDDecoder}
    z_window::Vector{Toric3DVertexMatchingDecoder}
    x_prefix::Vector{Toric3DVertexBPOSDDecoder}
    z_prefix::Vector{Toric3DVertexMatchingDecoder}
    final_x::Toric3DVertexBPOSDDecoder
    final_z::Toric3DVertexMatchingDecoder
    boundary::Toric3DVertexBoundaryMap
end
function _vertex_bposd(w,p,settings,ldpc,np,cache)
    prior=clamp(Float64(p),1e-12,1-1e-12)
    key=(size(w.matrix),String(UInt8.(vec(w.matrix))))
    decoder=get!(cache,key) do
        ldpc.BpOsdDecoder(_toric3d_numpy_matrix(np,w.matrix);
            error_rate=prior,max_iter=settings.max_iter,bp_method=settings.bp_method,
            ms_scaling_factor=settings.ms_scaling_factor,schedule=settings.schedule,
            omp_thread_count=settings.omp_thread_count,osd_method=settings.osd_method,osd_order=settings.osd_order)
    end
    Toric3DVertexBPOSDDecoder(w,Float64(p),prior,settings,decoder,np,"2.4.1")
end

"""Build every intermediate and final decoder; inputs are geometry and model parameters only."""
function build_toric3d_vertex_decoders(c::Toric3DVertexProjectorCode;p::Real,
        settings::Toric3DDecoderSettings=Toric3DDecoderSettings())
    isfinite(p) && 0<=p<=1 || throw(ArgumentError("p must lie in [0,1]"))
    ldpc=_toric3d_load_ldpc();version=_toric3d_module_version(ldpc,"ldpc")
    version=="2.4.1" || error("H+5 BP+OSD requires ldpc 2.4.1; found $version")
    np=PythonCall.pyimport("numpy");cache=Dict{Tuple{Tuple{Int,Int},String},PythonCall.Py}()
    bx(w)=_vertex_bposd(w,p,settings,ldpc,np,cache)
    xo=[build_toric3d_vertex_x_one_step_decoder(c,t) for t in 1:c.T-1]
    zo=[build_toric3d_vertex_z_one_step_decoder(c,t) for t in 1:c.T-1]
    xw=[bx(toric3d_vertex_x_complete_window(c,t)) for t in 1:c.T-1]
    zw=[_vertex_matching_decoder(toric3d_vertex_z_complete_window(c,t),:open) for t in 1:c.T-1]
    xp=[bx(toric3d_vertex_x_prefix_window(c,t)) for t in 1:c.T-1]
    zp=[_vertex_matching_decoder(toric3d_vertex_z_prefix_window(c,t),:open) for t in 1:c.T-1]
    qubits=collect(1:data_qubit_count(c))
    finalx=bx(_vertex_window(c,:x_only,:final,c.T,plaquette_checks(c),qubits,Int[],Int[]))
    finalz=_vertex_matching_decoder(_vertex_window(c,:z_only,:final,c.T,star_checks(c),qubits,Int[],Int[]),:closed)
    Toric3DVertexDecoderFamily(c,Float64(p),settings,xo,zo,xw,zw,xp,zp,finalx,finalz,toric3d_vertex_boundary_map(c))
end

function decode_toric3d_vertex_bposd(d::Toric3DVertexBPOSDDecoder,s::AbstractVector{Bool})::BitVector
    w=d.window
    length(s)==size(w.matrix,1) || throw(DimensionMismatch("plaquette syndrome length differs"))
    n=length(w.qubits)
    !any(s) && return falses(n)
    recovery=_toric3d_py_bits(d.numpy,d.bposd.decode(_toric3d_numpy_vector(d.numpy,s)),n)
    _vertex_validate_recovery(w,s,recovery)
end

"""Complete physical MWPM recovery, including its deterministic pair/path provenance."""
function decode_toric3d_vertex_z_complete(d::Toric3DVertexMatchingDecoder,s::AbstractVector{Bool})
    d.window.channel==:z_only && d.window.policy!=:one_step || throw(ArgumentError("not a complete Z decoder"))
    _vertex_matching_decision(d,s)
end

function toric3d_vertex_decoder_metadata(f::Toric3DVertexDecoderFamily)
    (p=f.p,p_eff=f.final_x.p_eff,settings=_toric3d_settings_metadata(f.settings),
     pymatching=f.final_z.version,ldpc=f.final_x.version,mwpm_weights="unit physical length",
     matching_tie_break="lexicographic sorted (a,b,physical-edge path) after minimum weight",
     neglected_correlations="independent physical-p priors omit gate-output propagation and circuit correlations")
end
