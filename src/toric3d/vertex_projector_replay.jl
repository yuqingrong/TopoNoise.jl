const TORIC3D_VERTEX_PROTOCOLS=(:none,:one_step,:window_complete,:full_prefix)
const TORIC3D_VERTEX_REPLAY_CHANNELS=(:x_only,:z_only,:both)

struct Toric3DVertexFeedbackHistory
    checkpoint::Int
    channel::Symbol
    window::Toric3DVertexFeedbackWindow
    syndrome::BitVector
    full_recovery::BitVector
    route::BitVector
    post_syndrome::BitVector
    paths::Vector{Toric3DVertexMatchedPath}
    signature::Vector{Tuple{Int,Int,Tuple}}
end
struct Toric3DVertexPreFinalState
    protocol::Symbol
    channel::Symbol
    before_boundary::PauliFrame
    boundary::Toric3DVertexBoundaryCompletion
    frame::PauliFrame
    history::Vector{Toric3DVertexFeedbackHistory}
    upper_parity_odd::Bool
end
struct Toric3DVertexReplayOutcome
    protocol::Symbol
    channel::Symbol
    before_boundary::PauliFrame
    boundary::Toric3DVertexBoundaryCompletion
    upper_parity_odd::Bool
    residual::PauliFrame
    history::Vector{Toric3DVertexFeedbackHistory}
    x_syndrome::BitVector
    z_syndrome::BitVector
    x_recovery::BitVector
    z_decision::Toric3DVertexMatchingDecision
    logical_x::NTuple{2,Bool}
    logical_z::NTuple{2,Bool}
    logical_x_failure::Bool
    logical_z_failure::Bool
    pauli_class::UInt8
    state_failure::Union{Nothing,Bool}
end

# Both sectors are read before either correction, including :both calls.
function _vertex_checkpoint_decisions(frame,f,protocol,channel,t)
    suffix=protocol==:one_step ? "one_step" : protocol==:window_complete ? "window" : "prefix"
    channels=channel==:both ? (:x_only,:z_only) : (channel,)
    decisions=Tuple{Toric3DVertexFeedbackWindow,BitVector,Toric3DVertexMatchingDecision}[]
    for ch in channels
        d=getproperty(f,Symbol(ch==:x_only ? "x_" : "z_",suffix))[t]
        w=d.window;bits=ch==:x_only ? frame.x : frame.z
        syndrome=_vertex_syndrome(w.checks,bits)
        if protocol==:one_step
            decision=ch==:x_only ? decode_toric3d_vertex_x_one_step(d,syndrome) : decode_toric3d_vertex_z_one_step(d,syndrome)
        elseif ch==:z_only
            decision=decode_toric3d_vertex_z_complete(d,syndrome)
        else
            recovery=decode_toric3d_vertex_bposd(d,syndrome)
            decision=Toric3DVertexMatchingDecision(Toric3DVertexMatchedPath[],recovery,copy(recovery),
                falses(length(syndrome)),Tuple{Int,Int,Tuple}[],"ldpc",d.version)
        end
        push!(decisions,(w,syndrome,decision))
    end
    decisions
end

function _vertex_final_recovery(frame::PauliFrame,f::Toric3DVertexDecoderFamily)
    xs=_vertex_syndrome(f.final_x.window.checks,frame.x)
    zs=_vertex_syndrome(f.final_z.window.checks,frame.z)
    x=decode_toric3d_vertex_bposd(f.final_x,xs)
    z=decode_toric3d_vertex_z_complete(f.final_z,zs)
    residual=PauliFrame(frame.x .⊻ x,frame.z .⊻ z.full_recovery)
    !any(_vertex_syndrome(f.final_x.window.checks,residual.x)) || error("final X recovery left syndrome")
    !any(_vertex_syndrome(f.final_z.window.checks,residual.z)) || error("final Z recovery left syndrome")
    lx=Tuple(_vertex_syndrome(toric3d_encoded_logical_z(f.code),residual.x))
    lz=Tuple(_vertex_syndrome(toric3d_encoded_logical_x(f.code),residual.z))
    cls=UInt8(sum(Int(b)<<(i-1) for (i,b) in enumerate((lx...,lz...))))
    (residual=residual,x_syndrome=xs,z_syndrome=zs,x_recovery=x,z_decision=z,
     logical_x=lx,logical_z=lz,logical_x_failure=any(lx),logical_z_failure=any(lz),pauli_class=cls)
end

"""Replay through upper-boundary completion, stopping before final recovery."""
function replay_toric3d_vertex_projector_pre_final(e::Toric3DVertexProjectorEncoder,
        r::Toric3DVertexProjectorFaultRecord,f::Toric3DVertexDecoderFamily;
        protocol::Symbol,channel::Symbol,record_history::Bool=true,_context=Dict{String,Any}())
    protocol in TORIC3D_VERTEX_PROTOCOLS || throw(ArgumentError("invalid H+5 feedback protocol"))
    channel in TORIC3D_VERTEX_REPLAY_CHANNELS || throw(ArgumentError("channel must be x_only, z_only, or both"))
    channel==:x_only && any(r.z) && throw(ArgumentError("X-only replay contains Z faults"))
    channel==:z_only && any(r.x) && throw(ArgumentError("Z-only replay contains X faults"))
    _vertex_validate_record(e,r)
    _vertex_dims(e.code)==_vertex_dims(f.code) || throw(DimensionMismatch("decoder and encoder dimensions differ"))
    c=e.code;m=c.Lu*c.Lv;n=data_qubit_count(c)
    frame=PauliFrame(falses(n),falses(n));history=Toric3DVertexFeedbackHistory[];first=1
    for t in 1:c.T
        propagate_toric3d_vertex_projector_segment!(frame,e,r,first,e.checkpoints[t]);first=e.checkpoints[t]+1
        (t==c.T || protocol==:none) && continue
        _context["frame"]=PauliFrame(copy(frame.x),copy(frame.z));_context["boundary_phase"]="pending";_context["checkpoint"]=t
        decisions=_vertex_checkpoint_decisions(frame,f,protocol,channel,t)
        for (w,s,decision) in decisions
            bits=w.channel==:x_only ? frame.x : frame.z
            all(q->data_qubit_coordinate(c,q)[2]<=t,w.qubits) || error("feedback correction accesses future data")
            _vertex_validate_recovery(w,s,decision.full_recovery)
            bits[w.qubits].⊻=decision.route
            post=_vertex_syndrome(w.checks,bits)
            post==decision.residual_syndrome || _vertex_invalid(w,s,decision.full_recovery;message="partial-route prediction disagrees with physical frame",route=decision.route,paths=[p.edges for p in decision.paths],observed_post_syndrome=post)
            if protocol!=:one_step
                checks=w.channel==:x_only ? plaquette_checks(c)[1:3m*t] : star_checks(c)[1:(t+1)*m]
                !any(_vertex_syndrome(checks,bits)) || _vertex_invalid(w,s,decision.full_recovery;message="complete feedback violated prepared-prefix checks",route=decision.route,paths=[p.edges for p in decision.paths],observed_post_syndrome=post)
            elseif w.channel==:z_only
                !any(post[w.retired_rows]) || _vertex_invalid(w,s,decision.full_recovery;message="one-step Z feedback left retired syndrome",route=decision.route,paths=[p.edges for p in decision.paths],observed_post_syndrome=post)
            end
            record_history && push!(history,Toric3DVertexFeedbackHistory(t,w.channel,w,s,decision.full_recovery,
                decision.route,post,decision.paths,decision.signature))
        end
    end
    before=PauliFrame(frame.x,frame.z)
    boundary=complete_toric3d_vertex_boundary(f.boundary,frame.z,r.upper_branch)
    _context["frame"]=PauliFrame(copy(frame.x),copy(frame.z));_context["boundary_phase"]="completed";_context["checkpoint"]=c.T
    frame.z.=boundary.z
    Toric3DVertexPreFinalState(protocol,channel,before,boundary,PauliFrame(frame.x,frame.z),history,
        isodd(count(boundary.observed)))
end

"""Replay one immutable shot, applying data corrections before subsequent gates."""
function _vertex_replay(e::Toric3DVertexProjectorEncoder,
        r::Toric3DVertexProjectorFaultRecord,f::Toric3DVertexDecoderFamily;
        protocol::Symbol,channel::Symbol,record_history::Bool=true,_context=Dict{String,Any}())
    pre=replay_toric3d_vertex_projector_pre_final(e,r,f;
        protocol,channel,record_history,_context)
    result=_vertex_final_recovery(pre.frame,f)
    Toric3DVertexReplayOutcome(protocol,channel,pre.before_boundary,pre.boundary,pre.upper_parity_odd,
        result.residual,pre.history,
        result.x_syndrome,result.z_syndrome,result.x_recovery,result.z_decision,result.logical_x,result.logical_z,
        result.logical_x_failure,result.logical_z_failure,result.pauli_class,
        channel==:x_only ? result.logical_x_failure : nothing)
end

function replay_toric3d_vertex_projector(e::Toric3DVertexProjectorEncoder,r::Toric3DVertexProjectorFaultRecord,f::Toric3DVertexDecoderFamily;kw...)
    context=Dict{String,Any}()
    try
        _vertex_replay(e,r,f;kw...,_context=context)
    catch err
        if err isa Toric3DVertexInvalidRecovery
            frame=pop!(context,"frame",PauliFrame(falses(data_qubit_count(e.code)),falses(data_qubit_count(e.code))))
            throw(Toric3DVertexInvalidRecovery(err.message,merge(context,_vertex_boundary_context(e,r,frame),err.context)))
        end
        rethrow()
    end
end
