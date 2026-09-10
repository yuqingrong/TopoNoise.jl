module TopoNoise

using ITensors: ITensors, ITensor, Index
using ITensorNetworks: ITensorNetwork
using LinearAlgebra: I, nullspace
using Random
using Statistics
import SparseArrays
import PythonCall
using Yao: H, X, Z, ResetTo, apply!, chain, control, matblock, measure!, put,
           repeat, subroutine, zero_state
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
       logical_x_support, logical_z_support,
       PlaquetteEncoder, rotated_planar_encoder, plaquette_blocks,
       gate_layers, yao_encoder, verify_encoder_tableau,
       CircuitPauliNoise, CircuitFaultStep, CircuitFaultRecord,
       sample_fault_record, with_pauli_fault, PauliFrame, SyndromeRecord,
       propagate_pauli_frame, measure_syndrome, syndrome_bits,
       sample_yao_syndrome,
       MatchingDecoders, build_matching_decoders, decode_logical_parities,
       LogicalFailurePoint, estimate_logical_failure,
       LogicalFailureScan, scan_logical_failure, write_logical_failure_csv,
       plot_logical_failure_scan, save_logical_failure_scan,
       ChannelFailureScan, ConstructionChannelComparison,
       scan_channel_logical_failure, run_construction_channel_comparison,
       comparison_series, channel_failure_count, channel_failure_rate,
       channel_failure_standard_error

include("toric_code_peps.jl")
include("local_gate.jl")
include("network_graph.jl")
include("sequential_circuit.jl")
include("visualization.jl")
include("rotated_planar/geometry.jl")
include("rotated_planar/encoder.jl")
include("rotated_planar/noise.jl")
include("rotated_planar/decoder.jl")
include("rotated_planar/scan.jl")
include("rotated_planar/comparison.jl")

end
