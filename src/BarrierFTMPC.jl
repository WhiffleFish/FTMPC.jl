module BarrierFTMPC

using LinearAlgebra
using StaticArrays
using ControlSystems
using BlockArrays
using SparseArrays
using OSQP

include("constants.jl")
include("linear.jl")
export LinearHexModel
export flip_z, trans_states

include("failure_linearization.jl")
export hover_control

include("nonlinear.jl")
export NonlinearHexModel

include("joint_dynamics.jl")
export joint_dynamics, make_joint, batch_dynamics, unbatch_states
export HexBatchDynamics

include("osqp.jl")
export OSQPFormulator, OSQPModel

# extras
include("extras.jl")
export basis, blkdiag

end
