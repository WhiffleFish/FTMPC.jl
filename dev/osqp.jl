using ControlSystems
using Plots
using LinearAlgebra
using OSQP
using BarrierFTMPC
using BlockArrays
const MPC = BarrierFTMPC

failures = [1]
T = 50
Δt = 0.1
sys = MPC.HexBatchDynamics(;failures, T, Δt)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[3] = -5

f = OSQPFormulator(sys;x_ref,Q=I(6)*1e-3)
model = OSQPModel(f, x0)
res = OSQP.solve!(model)

x, u = res.x[1:T*12], res.x[T*12+1:end]

X = MPC.unbatch_states(x, 12)
p1 = plot(trans_states(flip_z(X))', labels=reshape(MPC.STATE_LABELS[MPC.TRANSLATIONAL_STATES], 1,6))
p2 = plot(reshape(u, 6, length(u) ÷ 6)', labels=reshape(["u$i" for i ∈ 1:6],1,6))
plot(p1, p2)


m1, m2 = MPC.unbatch_and_disjoint(res.x[1:12*nm*T], nm, T)
m1

6*(T-1)*nm
u = res.x[12*nm*T+1:end]
