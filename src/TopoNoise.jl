module TopoNoise

using ITensors: ITensors, ITensor, Index
using ITensorNetworks: ITensorNetwork
using LinearAlgebra: I, nullspace
import Random
import Statistics
import Yao
using Yao: H, chain, matblock, repeat, subroutine
import CairoMakie
import ITensorNetworks

export ToricCodePEPS, physicalinds, toric_code_local_tensor, toric_code_peps,
       toric_code_local_gate, unitary_completion, peps_graph,
       ToricCodeSequentialCircuit, sequential_circuit_graph, circuit_layers,
       yao_unitary, yao_circuit,
       log2_postselection_probability, postselection_probability,
       ToricCodeTrajectoryModel, VirtualBondErrors, ToricCodeTrajectory,
       sample_virtual_errors, sample_trajectory, sample_yao_trajectory,
       TrajectoryObservables, bond_mismatches, trajectory_observables,
       TrajectoryScanPoint, TrajectoryScan, CriticalCrossing,
       scan_trajectories, estimate_crossings,
       plot_peps_graph, plot_sequential_circuit, plot_trajectory_scan

include("toric_code_peps.jl")
include("local_gate.jl")
include("network_graph.jl")
include("sequential_circuit.jl")
include("trajectory.jl")
include("trajectory_analysis.jl")
include("visualization.jl")

end
