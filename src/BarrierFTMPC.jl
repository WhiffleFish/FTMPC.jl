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

include(joinpath("models", "models.jl"))

include("joint_dynamics.jl")
export joint_dynamics, make_joint, batch_dynamics, unbatch_states

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

include("envs.jl")

# extras
include("extras.jl")
export basis, blkdiag

end
