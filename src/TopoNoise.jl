module TopoNoise

using ITensors: ITensors, ITensor, Index
using ITensorNetworks: ITensorNetwork
using LinearAlgebra: I, nullspace, rank
import Random
import Statistics
import Optim
import Yao
using Yao: H, chain, matblock, repeat, subroutine
import CairoMakie
import ITensorNetworks
using PythonCall: Py, pyimport, pyconvert, pybuiltins, pylist

export ToricCodePEPS, physicalinds, toric_code_local_tensor, toric_code_peps,
       toric_code_local_gate, unitary_completion, peps_graph,
       ToricCodeSequentialCircuit, sequential_circuit_graph, circuit_layers,
       yao_unitary, yao_circuit,
       log2_postselection_probability, postselection_probability,
       ToricCodeTrajectoryModel, VirtualBondErrors, ToricCodeTrajectory,
       sample_virtual_errors, sample_trajectory, sample_yao_trajectory,
       TrajectoryObservables, bond_mismatches, trajectory_observables,
       MarginalSpinObservables, marginal_spin_observables,
       TrajectoryScanPoint, MarginalSpinScanPoint, TrajectoryScan,
       CriticalCrossing, scan_trajectories, estimate_crossings,
       estimate_binder_crossings,
       RotatedCodeCapacityModel, RotatedCodeCapacityScanPoint,
       RotatedCodeCapacityScan, rotated_code_capacity_circuit,
       estimate_rotated_code_capacity, scan_rotated_code_capacity,
       estimate_rotated_code_crossings, plot_rotated_code_capacity,
       IsometricPlanarCodeModel, IsometricPlanarVirtualErrors,
       IsometricPlanarSyndrome, sample_isometric_planar_virtual_errors,
       isometric_planar_syndrome, decode_isometric_planar_syndrome,
       isometric_planar_logical_failure, IsometricPlanarEncoder,
       IsometricPlanarCheckLayer, isometric_planar_encoder,
       isometric_planar_check_layer,
       IsometricPlanarYaoTrajectory, sample_isometric_planar_yao_trajectory,
       isometric_planar_bond_mismatches,
       IsometricPlanarCapacityPoint, IsometricPlanarCapacityScan,
       estimate_isometric_planar_capacity, scan_isometric_planar_capacity,
       estimate_isometric_planar_crossings, IsometricPlanarScalingFit,
       fit_isometric_planar_scaling, IsometricPlanarScalingDiagnostics,
       diagnose_isometric_planar_scaling,
       plot_isometric_planar_capacity,
       plot_isometric_planar_scaling_diagnostics,
       Correction, DecodedTrajectory, plaquette_syndrome, decode_uf,
       logical_failure, decode_trajectory,
       plot_peps_graph, plot_sequential_circuit, plot_trajectory_scan

include("toric_code_peps.jl")
include("local_gate.jl")
include("network_graph.jl")
include("sequential_circuit.jl")
include("trajectory.jl")
include("code_capacity.jl")
include("decoder.jl")
include("trajectory_analysis.jl")
include("isometric_planar_code.jl")
include("visualization.jl")

end
