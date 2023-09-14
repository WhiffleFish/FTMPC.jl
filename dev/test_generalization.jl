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

models = [MPC.HexCTLinearModel(failure) for failure in failures]
sys = MPC.BatchDynamics(models; T, Δt, u_bounds)

x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= [-0.4, 0.4, 5]

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
a, info = MPC.action_info(unit_planner, rand(12)*1e-2)

unit_sim = Simulator(MPC.HexIMM(;Δt), unit_planner, x0=x0, T=100, failure=MPC.FixedFailure(10,3;instant_update=true))
unit_hist = simulate(unit_sim)
plot(unit_hist, lw=2)
