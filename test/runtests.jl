using ControlSystems
using Test
using LinearAlgebra
using BarrierFTMPC
const MPC = BarrierFTMPC

@testset "lqr" begin
    include("lqr.jl")
end

@testset "batch" begin
    include("batch.jl")
end
