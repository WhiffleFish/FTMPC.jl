using ControlSystems
using Plots
using LinearAlgebra
using OSQP
using BarrierFTMPC
const MPC = BarrierFTMPC

failures = [0,1,2,3]
T = 100
Δt = 0.1
# u_bounds = (0.,10.)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5


f = OSQPFormulator(sys;x_ref,Q=I(6)*1e-2, polish=true)
model = OSQPModel(f, x0)
res = OSQP.solve!(model)

_res = MPC.HexOSQPResults(f, res)

p1 = plot(_res.t,
    trans_states(flip_z(_res.X[1]))',
    labels=permutedims(MPC.STATE_LABELS[MPC.TRANSLATIONAL_STATES]),
    lw=2
)
p2 = plot(_res.t[1:end-1],_res.U[1]', lw=2, labels=permutedims(["u$i" for i ∈ 1:6]))

plot(p1,p2)

plot(_res.t,
    trans_states(flip_z(_res.X[1]))',
    labels=permutedims(MPC.STATE_LABELS[MPC.TRANSLATIONAL_STATES]),
    lw=2
)
plot(_res.t,
    trans_states(flip_z(_res.X[2]))',
    labels=permutedims(MPC.STATE_LABELS[MPC.TRANSLATIONAL_STATES]),
    lw=2
)

plot(_res.t[1:end-1],_res.U[1]', lw=2, labels=permutedims(["u$i" for i ∈ 1:6]))
plot(_res.t[1:end-1],_res.U[2]', lw=2)
