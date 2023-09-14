using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Plots

default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")
constraints = MPC.ElevatorShaft(h=6, γ=1e-1)
failures = 0:6
T = 10
Δt = 0.05

u_bounds = (-Inf,Inf)
nm = length(failures)

# FIXME: should be different stuff here
models = [MPC.HexCTLinearModel(failure) for failure in failures]
sys = MPC.BatchDynamics(models; T, Δt, u_bounds)

x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= [-0.4, 0.4, 5]

q_vec = zeros(12)
q_vec[1:3] .= 1.
Q = diagm(q_vec)

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    Q=I(12),
    R=I(6)*1e-3,
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-5,
    eps_rel = 1e-5,
    verbose = false,
    max_iter= 50_000
)
model = JuMPModel(f, x0)
unit_planner = MPC.FTMPCPlanner(model, f, 1)

using BenchmarkTools
@btime MPC.action(unit_planner, rand(12)*0.01)
