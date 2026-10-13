"""Diagonal H+5 slab: open transfer depth T, periodic transverse periods Lu,Lv."""
struct Toric3DVertexProjectorCode
    T::Int
    Lu::Int
    Lv::Int
    encoded::Toric3DEncodedSlice
    function Toric3DVertexProjectorCode(T::Integer, Lu::Integer, Lv::Integer)
        encoded = Toric3DEncodedSlice(T, Lu, Lv)
        new(Int(T), Int(Lu), Int(Lv), encoded)
    end
end

data_qubit_count(c::Toric3DVertexProjectorCode) = data_qubit_count(c.encoded)
data_qubit_index(c::Toric3DVertexProjectorCode, args...) = data_qubit_index(c.encoded,args...)
data_qubit_coordinate(c::Toric3DVertexProjectorCode,q::Integer) = data_qubit_coordinate(c.encoded,q)
data_qubit_coordinates(c::Toric3DVertexProjectorCode) = data_qubit_coordinates(c.encoded)
star_checks(c::Toric3DVertexProjectorCode) = star_checks(c.encoded)
plaquette_checks(c::Toric3DVertexProjectorCode) = plaquette_checks(c.encoded)
star_check_matrix(c::Toric3DVertexProjectorCode; sparse::Bool=false) = star_check_matrix(c.encoded;sparse)
plaquette_check_matrix(c::Toric3DVertexProjectorCode; sparse::Bool=false) = plaquette_check_matrix(c.encoded;sparse)
toric3d_encoded_logical_x(c::Toric3DVertexProjectorCode) = toric3d_encoded_logical_x(c.encoded)
toric3d_encoded_logical_z(c::Toric3DVertexProjectorCode) = toric3d_encoded_logical_z(c.encoded)
toric3d_slice_qubits(c::Toric3DVertexProjectorCode,t::Integer) = toric3d_slice_qubits(c.encoded.lattice,t)
toric3d_slice_checks(c::Toric3DVertexProjectorCode,t::Integer) = toric3d_slice_checks(c.encoded.lattice,t)

_vertex_dims(c::Toric3DVertexProjectorCode) = (c.T,c.Lu,c.Lv)
_vertex_syndrome(h::AbstractMatrix{Bool},v::AbstractVector{Bool}) =
    BitVector(isodd(count(@view(h[r,:]) .& v)) for r in axes(h,1))
_vertex_syndrome(checks::AbstractVector,v::AbstractVector{Bool}) =
    BitVector(isodd(count(v[s])) for s in checks)
