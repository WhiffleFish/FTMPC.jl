using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using COSMO
using LinearAlgebra
using Plots
default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

failures = [0]
T = 50
Δt = 0.1
u_bounds = (0.1,15.)
# u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5
ss = LinearHexModel(0)

f = BarrierJuMPFormulator(sys, OSQP.Optimizer;x_ref,Q=I(6)*1e-2,verbose=false)
model = JuMPModel(f, x0)

planner = FTMPCPlanner(model, f)
sim = MPC.Simulator(LinearHexModel(0), planner, x0=zeros(12), T=50)
hist = simulate(sim)

plot(hist, lw=2)

## with barriers

constraints = [
    LinearConstraint(basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(basis(12, 2)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 2)*1, 1, 1e-1)
]

failures = [0,1]
T = 10
Δt = 0.1

u_bounds = (0.0,15.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -10

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    # Q=(I(6)*1e-1,1),
    # R=(I(6)*1e-1,1),
    Q=I(6)*1e-1,
    R=I(6)*1e-1,
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-4,
    eps_rel = 1e-4,
    verbose = true,
    max_iter= 50_000
)
model = JuMPModel(f, x0)
# MPC.set_consensus_horizon(model, f, 1)
planner = FTMPCPlanner(model, f)
sim = Simulator(LinearHexModel(0), planner, x0=x0, T=200)
hist = simulate(sim)

plot(hist, lw=2)
