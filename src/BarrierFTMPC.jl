module BarrierFTMPC

using LinearAlgebra
using StaticArrays
using ControlSystemsBase
using BlockArrays
using SparseArrays
using ProgressMeter
using RecipesBase
using Random
using Distributions: pdf, MvNormal
using JuMP
using Plots

include("constants.jl")
include("linear.jl")
export LinearHexModel
export flip_z, trans_states, pos_states

include("failure_linearization.jl")
export hover_control

include("joint_dynamics.jl")
export joint_dynamics, make_joint, batch_dynamics, unbatch_states
export HexBatchDynamics

include("barriers.jl")
export BarrierJuMPFormulator, LinearConstraint, JuMPModel

include("binary_search.jl")
export binary_search_max

include("planner.jl")
export FTMPCPlanner

include("imm.jl")
export IMM

include("simulator.jl")
export Simulator, simulate

include("plotting.jl")

# extras
include("extras.jl")
export basis, blkdiag

include("envs.jl")

end
