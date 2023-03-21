module BarrierFTMPC

using LinearAlgebra
using StaticArrays
using ControlSystemsBase
using BlockArrays
using SparseArrays
using ProgressMeter
using RecipesBase
using JuMP
using OSQP

include("constants.jl")
include("linear.jl")
export LinearHexModel
export flip_z, trans_states

include("binary_search.jl")
export binary_search_max

include("failure_linearization.jl")
export hover_control

include("joint_dynamics.jl")
export joint_dynamics, make_joint, batch_dynamics, unbatch_states
export HexBatchDynamics

include("osqp.jl")
export OSQPFormulator, OSQPModel

include("jump.jl")
export JuMPFormulator, JuMPModel

include("barriers.jl")
export BarrierJuMPFormulator, LinearConstraint

include("planner.jl")
export FTMPCPlanner

include("simulator.jl")
export Simulator, simulate

# extras
include("extras.jl")
export basis, blkdiag

end
