using ControlSystems
using Test
using LinearAlgebra
using BarrierFTMPC
using OSQP
const MPC = BarrierFTMPC

@testset "lqr" begin
    include("lqr.jl")
end

@testset "batch" begin
    include("batch.jl")
end

@testset "osqp" begin
    include("osqp.jl")
end
