using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using OSQP
using LinearAlgebra
using Plots

constraints = [
    LinearConstraint(basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 3)*1, 1, 1e-1)
]

failures = 0:6
T = 20
Δt = 0.1
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= 5

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    P = (I(12),1),#(I(12),1),
    Q = (I(6)*1e-2,1),# (I(6)*1e-2,1),
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-5,
    eps_rel = 1e-5,
    verbose = false
)

model = JuMPModel(f, x0)
planner = MPC.ConsensusSearchPlanner(model, f)
sim = Simulator(LinearHexModel(0), planner, x0=x0, T=50)
hist = simulate(sim)
plot(hist)
