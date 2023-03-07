module BarrierFTMPC

using LinearAlgebra
using StaticArrays
using ControlSystems
using BlockArrays
# using OSQP

include("constants.jl")
include("linear.jl")
export LinearHexModel
export flip_z, trans_states

include("failure_linearization.jl")
export hover_control

include("nonlinear.jl")
export NonlinearHexModel

include("belief.jl")

# extras
export basis

end
