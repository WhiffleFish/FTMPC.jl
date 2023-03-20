using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Plots

failures = [0,1]
T = 50
Δt = 0.1
u_bounds = (0.,12.)
# u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5
ss = LinearHexModel(0)

f = JuMPFormulator(sys, OSQP.Optimizer;x_ref,Q=I(6)*1e-2, polish=false, verbose=false)
model = JuMPModel(f, x0)
optimize!(model) # warm start

planner = FTMPCPlanner(model, f)
sim = Simulator(LinearHexModel(0), planner, x0=zeros(12), T=50)
hist = simulate(sim)

plot(hist, lw=2)

## with barriers

A_constraint = zeros(12)
A_constraint[2] = -1
b_constraint = 1
γ_constraint = 5e-2

failures = [0,1]
T = 100
Δt = 0.1

# u_bounds = (0.,10.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[2] = 10

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    Q=I(6)*1e-2,
    A_constraint,
    b_constraint,
    γ_constraint,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-4,
    eps_rel = 1e-4,
    verbose = false,
    time_limit = Δt
)
model = JuMPModel(f, x0)
optimize!(model)

planner = FTMPCPlanner(model, f)
sim = Simulator(LinearHexModel(0), planner, x0=zeros(12), T=100)
hist = simulate(sim)

plot(hist, lw=2)
