module TopoNoise

using ITensors: ITensors, ITensor, Index
using ITensorNetworks: ITensorNetwork
using Distributions: Binomial
using LinearAlgebra: I, nullspace
using Random
using Statistics
import TOML
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
       channel_failure_standard_error,
       PairCrossing, ThresholdFit, FittedConstructionChannelComparison,
       fit_channel_threshold, fit_construction_channel_comparison,
       scaled_error_rate,
       write_comparison_raw_csv, write_comparison_fit_csv,
       plot_construction_channel_comparison, plot_construction_channel_panel,
       save_construction_channel_comparison

export ToricCode3D, star_checks, plaquette_checks, star_check_matrix, plaquette_check_matrix,
       toric3d_prefix, Toric3DEncoder, toric3d_encoder, verify_toric3d_encoder,
       Toric3DFaultRecord, sample_toric3d_faults, with_toric3d_fault,
       propagate_toric3d_segment!, toric3d_syndrome, toric3d_logical_parities,
       Toric3DDecoderSettings, build_toric3d_decoder, decode_toric3d,
       toric3d_decoder_metadata, Toric3DMetricEstimate, Toric3DComparisonPoint,
       Toric3DComparison, replay_toric3d_feedback, compare_toric3d_feedback,
       toric3d_z_lower_bound, plot_toric3d_comparison, save_toric3d_comparison

export Toric3DSliceLattice, Toric3DSliceEncoder, toric3d_slice_qubits,
       toric3d_slice_boundary_checks, toric3d_slice_checks, toric3d_slice_check_matrix,
       toric3d_slice_encoder, verify_toric3d_slice_encoder,
       build_toric3d_slice_decoder, decode_toric3d_slice, replay_toric3d_slice,
       Toric3DSliceMetricEstimate, Toric3DSliceComparisonPoint, Toric3DSliceComparison,
       toric3d_slice_diagnostics, compare_toric3d_slice_feedback,
       plot_toric3d_slice_comparison, save_toric3d_slice_comparison

export Toric3DEncodedSlice, toric3d_encoded_logical_x, toric3d_encoded_logical_z,
       toric3d_encoded_seed_checks, toric3d_encoded_encoder, verify_toric3d_encoded_encoder,
       verify_toric3d_encoded_zero_encoder,
       build_toric3d_z_fault_model, sample_toric3d_z_frames,
       build_toric3d_z_threshold_decoder, decode_toric3d_z_batch,
       toric3d_z_threshold_counts, scan_toric3d_z_threshold, save_toric3d_z_threshold,
       Toric3DXThresholdDecoder, build_toric3d_x_threshold_decoder,
       decode_toric3d_x, toric3d_x_threshold_outcome,
       replay_toric3d_x_feedback, toric3d_x_feedback_pair,
       Toric3DXFeedbackPoint, Toric3DXFeedbackScan,
       compare_toric3d_x_feedback_threshold

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
include("rotated_planar/scaling.jl")
include("rotated_planar/comparison_plot.jl")

include("toric3d/geometry.jl")
include("toric3d/encoder.jl")
include("toric3d/noise.jl")
include("toric3d/decoder.jl")
include("toric3d/scan.jl")
include("toric3d/output.jl")
include("toric3d/slice_geometry.jl")
include("toric3d/slice_encoder.jl")
include("toric3d/slice_feedback.jl")
include("toric3d/slice_scan.jl")
include("toric3d/encoded_slice_geometry.jl")
include("toric3d/z_threshold.jl")
include("toric3d/x_feedback_threshold.jl")
include("toric3d/x_feedback_scan.jl")

end
