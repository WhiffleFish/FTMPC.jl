using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using OSQP
using LinearAlgebra
using Plots

failures = [0,1]
T = 100
Δt = 0.1
# u_bounds = (0.,10.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5
ss = LinearHexModel(0)

f = JuMPFormulator(sys, OSQP.Optimizer;x_ref,Q=I(6)*1e-2, polish=false, verbose=false)
model = JuMPModel(f, x0)

planner = FTMPCPlanner(model, f)
sim = Simulator(LinearHexModel(0), planner, zeros(12), 50)
hist = simulate(sim)

plot(hist.x')
plot(hist.u')
