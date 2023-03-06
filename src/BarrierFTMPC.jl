module BarrierFTMPC

using LinearAlgebra
using StaticArrays
using ControlSystems

include("constants.jl")
include("linear.jl")
export LinearHexModel
export flip_z, trans_states

include("nonlinear.jl")
export NonlinearHexModel

# extras
export basis

end
