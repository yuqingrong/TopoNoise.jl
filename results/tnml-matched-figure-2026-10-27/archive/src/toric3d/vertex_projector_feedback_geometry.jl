"""Measured full-code supports and the restricted data columns allowed to change."""
struct Toric3DVertexFeedbackWindow
    code::Toric3DVertexProjectorCode
    channel::Symbol
    policy::Symbol
    checkpoint::Int
    checks::Vector{Vector{Int}}
    qubits::Vector{Int}
    matrix::BitMatrix
    retired_rows::Vector{Int}
    frontier_rows::Vector{Int}
end
function _vertex_window(c,channel,policy,t,checks,qubits,retired,frontier)
    Toric3DVertexFeedbackWindow(c,channel,policy,Int(t),checks,qubits,
        _toric3d_restricted_check_matrix(checks,qubits),collect(retired),collect(frontier))
end
function _vertex_checkpoint(c,t)
    1<=t<c.T || throw(ArgumentError("intermediate feedback requires 1 <= t < T"))
end
function toric3d_vertex_x_slice_window(c::Toric3DVertexProjectorCode,t::Integer)
    _vertex_checkpoint(c,t)
    checks=toric3d_slice_checks(c,t)
    _vertex_window(c,:x_only,:one_step,t,checks,toric3d_slice_qubits(c,t),Int[],1:length(checks))
end
function toric3d_vertex_x_complete_window(c::Toric3DVertexProjectorCode,t::Integer)
    _vertex_checkpoint(c,t);m=c.Lu*c.Lv
    checks=plaquette_checks(c)[max(1,3m*(t-2)+1):3m*t]
    qubits=vcat(toric3d_slice_qubits(c,t-1),toric3d_slice_qubits(c,t))
    _vertex_window(c,:x_only,:window_complete,t,checks,qubits,t==1 ? Int[] : (1:3m),
        t==1 ? (1:3m) : (3m+1:6m))
end
function toric3d_vertex_x_prefix_window(c::Toric3DVertexProjectorCode,t::Integer)
    _vertex_checkpoint(c,t);m=c.Lu*c.Lv
    _vertex_window(c,:x_only,:full_prefix,t,plaquette_checks(c)[1:3m*t],
        collect(1:3m*(t+1)),1:3m*(t-1),3m*(t-1)+1:3m*t)
end
function _vertex_z_window(c,t,policy)
    _vertex_checkpoint(c,t);m=c.Lu*c.Lv
    checks=star_checks(c)[(t-1)*m+1:(t+1)*m]
    qubits=vcat(toric3d_slice_qubits(c,t-1),toric3d_slice_qubits(c,t))
    _vertex_window(c,:z_only,policy,t,checks,qubits,1:m,m+1:2m)
end
toric3d_vertex_z_strip_window(c::Toric3DVertexProjectorCode,t::Integer)=_vertex_z_window(c,t,:one_step)
toric3d_vertex_z_complete_window(c::Toric3DVertexProjectorCode,t::Integer)=_vertex_z_window(c,t,:window_complete)
function toric3d_vertex_z_prefix_window(c::Toric3DVertexProjectorCode,t::Integer)
    _vertex_checkpoint(c,t);m=c.Lu*c.Lv
    _vertex_window(c,:z_only,:full_prefix,t,star_checks(c)[1:(t+1)*m],
        collect(1:3m*(t+1)),1:t*m,t*m+1:(t+1)*m)
end
