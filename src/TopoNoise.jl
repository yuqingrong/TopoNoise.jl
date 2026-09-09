module TopoNoise

using ITensors: ITensors, ITensor, Index
using ITensorNetworks: ITensorNetwork
using LinearAlgebra: I, nullspace
import SparseArrays
using Yao: H, chain, matblock, repeat, subroutine
import CairoMakie
import ITensorNetworks

export ToricCodePEPS, physicalinds, toric_code_local_tensor, toric_code_peps,
       toric_code_local_gate, unitary_completion, peps_graph,
       ToricCodeSequentialCircuit, sequential_circuit_graph, circuit_layers,
       yao_unitary, yao_circuit,
       log2_postselection_probability, postselection_probability,
       plot_peps_graph, plot_sequential_circuit,
       RotatedPlanarCode, StabilizerCheck, distance, boundary_orientation,
       data_qubit_count, data_qubit_index, data_qubit_coordinate,
       data_qubit_coordinates, a_s_checks, b_p_checks, stabilizers,
       a_s_check_matrix, b_p_check_matrix, stabilizer_check_matrix,
       logical_x_support, logical_z_support

include("toric_code_peps.jl")
include("local_gate.jl")
include("network_graph.jl")
include("sequential_circuit.jl")
include("visualization.jl")
include("rotated_planar/geometry.jl")

end
