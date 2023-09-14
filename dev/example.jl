using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Plots
default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

constraints = [
    LinearConstraint(basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 3)*0.5, 1, 1e-1)
]

failures = [0,1]
T = 50
Δt = 0.1

sys = MPC.HexBatchDynamics(;failures, T, Δt)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5
ss = LinearHexModel(0)

f = BarrierJuMPFormulator(
    sys, 
    OSQP.Optimizer;
    x_ref,
    constraints,
    Q = I(6)*1e-2, 
    verbose=false
)
model = JuMPModel(f, x0)
optimize!(model)

res = MPC.OSQPResults(f, model)
p1 = plot(res.t,
    trans_states(flip_z(res.X[1]))',
    labels=permutedims(MPC.STATE_LABELS[MPC.TRANSLATIONAL_STATES]),
    lw=2
)
p2 = plot(
    res.t[1:end-1],
    res.U[1]', 
    lw=2, 
    labels=permutedims(["u$i" for i ∈ 1:6])
)
plot(p1,p2)
planner = FTMPCPlanner(model, f)
sim = Simulator(LinearHexModel(0), planner, x0=x0, T=100)
hist = simulate(sim)

plot(hist, lw=2)
