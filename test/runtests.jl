using ControlSystemsBase
using Test
using LinearAlgebra
using BarrierFTMPC
using OSQP
const MPC = BarrierFTMPC
using JuMP

@testset "batch" begin
    include("batch.jl")
end

@testset "barrier" begin
    include("barrier.jl")
end

@testset "simulate" begin
    include("simulate.jl")
end

@testset "b-search" begin
    include("binary_search.jl")
end
