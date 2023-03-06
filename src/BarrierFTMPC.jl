module BarrierFTMPC

using LinearAlgebra
using StaticArrays
using ControlSystems

include(joinpath(@__DIR__, "hexmodel.jl"))
export LinearHexModel

end
