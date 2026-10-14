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
       Toric3DEncodedPrefix, toric3d_encoded_prefix,
       Toric3DZFeedbackWindow, toric3d_z_feedback_window, toric3d_z_window_syndrome,
       build_toric3d_z_fault_model, sample_toric3d_z_frames,
       build_toric3d_z_threshold_decoder, decode_toric3d_z_batch,
       toric3d_z_threshold_outcome, toric3d_z_threshold_counts,
       Toric3DZFeedbackDecoder, Toric3DZFeedbackDecoderFamily,
       build_toric3d_z_feedback_decoders, decode_toric3d_z_feedback,
       decode_toric3d_z_feedback_batch,
       Toric3DDiagonalZFaultRecord, sample_toric3d_diagonal_z_faults,
       replay_toric3d_z_strength, toric3d_z_strength_triplet,
       Toric3DDiagonalZFaultBatch, sample_toric3d_diagonal_z_fault_batch,
       toric3d_z_strength_batch,
       scan_toric3d_z_threshold, save_toric3d_z_threshold,
       Toric3DXThresholdDecoder, build_toric3d_x_threshold_decoder,
       decode_toric3d_x, toric3d_x_threshold_outcome,
       replay_toric3d_x_feedback, toric3d_x_feedback_pair,
       Toric3DXFeedbackDecoderFamily, build_toric3d_x_feedback_decoders,
       replay_toric3d_x_strength, toric3d_x_strength_triplet,
       Toric3DProtocolEstimate, Toric3DPairEstimate, Toric3DTripleMetricEstimate,
       Toric3DDiagonalFeedbackPoint, Toric3DDiagonalFeedbackScan,
       compare_toric3d_diagonal_feedback, plot_toric3d_diagonal_feedback,
       save_toric3d_diagonal_feedback,
       Toric3DXFeedbackPoint, Toric3DXFeedbackScan,
       compare_toric3d_x_feedback_threshold,
       plot_toric3d_x_feedback_threshold, save_toric3d_x_feedback_threshold

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
include("toric3d/feedback_strength_geometry.jl")
include("toric3d/z_threshold.jl")
include("toric3d/feedback_strength_decoders.jl")
include("toric3d/feedback_strength_z.jl")
include("toric3d/feedback_strength_z_batch.jl")
include("toric3d/x_feedback_threshold.jl")
include("toric3d/feedback_strength_x.jl")
include("toric3d/feedback_strength_scan.jl")
include("toric3d/feedback_strength_output.jl")
include("toric3d/x_feedback_scan.jl")
include("toric3d/x_feedback_output.jl")

export Toric3DVertexProjectorCode, Toric3DVertexProjectorEncoder,
       toric3d_vertex_projector_encoder, verify_toric3d_vertex_projector_encoder
include("toric3d/vertex_projector_geometry.jl")
include("toric3d/vertex_projector_encoder.jl")

export Toric3DVertexBoundaryMap, toric3d_vertex_boundary_map,
       verify_toric3d_vertex_lower_boundary, toric3d_vertex_boundary_correction,
       sample_toric3d_vertex_upper_branch, Toric3DVertexBoundaryCompletion,
       complete_toric3d_vertex_boundary
include("toric3d/vertex_projector_boundary.jl")

export Toric3DVertexProjectorFaultRecord, sample_toric3d_vertex_projector_record,
       with_toric3d_vertex_projector_fault, propagate_toric3d_vertex_projector_segment!
include("toric3d/vertex_projector_noise.jl")

export Toric3DVertexFeedbackWindow, toric3d_vertex_x_slice_window,
       toric3d_vertex_x_complete_window, toric3d_vertex_x_prefix_window,
       toric3d_vertex_z_strip_window, toric3d_vertex_z_complete_window,
       toric3d_vertex_z_prefix_window
include("toric3d/vertex_projector_feedback_geometry.jl")

export Toric3DVertexMatchedPath, Toric3DVertexMatchingDecision,
       build_toric3d_vertex_x_one_step_decoder, build_toric3d_vertex_z_one_step_decoder,
       decode_toric3d_vertex_x_one_step, decode_toric3d_vertex_z_one_step
include("toric3d/vertex_projector_matching.jl")

export Toric3DVertexBPOSDDecoder, Toric3DVertexDecoderFamily,
       build_toric3d_vertex_decoders, decode_toric3d_vertex_bposd,
       decode_toric3d_vertex_z_complete, toric3d_vertex_decoder_metadata
include("toric3d/vertex_projector_decoders.jl")

export TORIC3D_VERTEX_PROTOCOLS, TORIC3D_VERTEX_REPLAY_CHANNELS,
       Toric3DVertexFeedbackHistory, Toric3DVertexPreFinalState, Toric3DVertexReplayOutcome,
       replay_toric3d_vertex_projector_pre_final, replay_toric3d_vertex_projector
include("toric3d/vertex_projector_replay.jl")

export Toric3DVertexZFaultSignature, Toric3DVertexZObservableGroup,
       Toric3DVertexZParitySummary, Toric3DVertexZAmbiguityAudit,
       toric3d_vertex_z_fault_signatures, audit_toric3d_vertex_z_boundary
include("toric3d/vertex_projector_z_audit.jl")

export Toric3DVertexCircuitZDecoder, Toric3DVertexCircuitZDecision,
       Toric3DVertexZControlOutcome, Toric3DVertexParityErasureEstimate,
       build_toric3d_vertex_circuit_z_decoder, decode_toric3d_vertex_circuit_z,
       finalize_toric3d_vertex_z_control, toric3d_vertex_parity_erasure_estimate
include("toric3d/vertex_projector_z_control.jl")

export Toric3DVertexScanGrid, Toric3DVertexProtocolEstimate, Toric3DVertexPairEstimate,
       Toric3DVertexQuadMetricEstimate, Toric3DVertexFeedbackPoint, Toric3DVertexFeedbackScan,
       compare_toric3d_vertex_projector_feedback, Toric3DVertexInvalidRecovery
include("toric3d/vertex_projector_scan.jl")

export write_toric3d_vertex_projector_feedback_csv, plot_toric3d_vertex_projector_feedback,
       save_toric3d_vertex_projector_feedback, Toric3DVertexProductionSelection,
       load_toric3d_vertex_production_selection
include("toric3d/vertex_projector_output.jl")

end
