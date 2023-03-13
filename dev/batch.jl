using ControlSystems
using Plots
using LinearAlgebra
using BarrierFTMPC
using BlockArrays
const MPC = BarrierFTMPC

failures = [0]
T = 10
Δt = 0.1
sys = MPC.HexBatchDynamics(;failures, T, Δt)
x0 = zeros(12)
sys.A
sys.A*x0
