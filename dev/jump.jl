using ControlSystems
using Plots
using LinearAlgebra
using OSQP
using BarrierFTMPC
using JuMP
const MPC = BarrierFTMPC

failures = [0,1]
T = 100
Δt = 0.1
u_bounds = (0.,10.)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5



f = JuMPFormulator(sys, OSQP.Optimizer;x_ref,Q=I(6)*1e-2, polish=true)
model = JuMPModel(f, x0)
optimize!(model)
res = MPC.HexOSQPResults(f, model)

p1 = plot(res.t,
    trans_states(flip_z(res.X[1]))',
    labels=permutedims(MPC.STATE_LABELS[MPC.TRANSLATIONAL_STATES]),
    lw=2
)
p2 = plot(res.t[1:end-1],res.U[1]', lw=2, labels=permutedims(["u$i" for i ∈ 1:6]))
plot(p1,p2)
