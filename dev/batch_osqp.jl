using ControlSystems
using Plots
using LinearAlgebra
using OSQP
using BarrierFTMPC
using BlockArrays
const MPC = BarrierFTMPC

failures = [0,1,2]
T = 10
Δt = 0.1
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[3] = -5


f = OSQPFormulator(sys;x_ref,Q=I(6)*1e-3)
model = OSQPModel(f, x0)
res = OSQP.solve!(model)

_res = MPC.HexOSQPResults(f, res)
plot(_res.X[1]')
plot(_res.X[2]')

plot(_res.U[1]')
plot(_res.U[2]')

x, u = res.x[1:T*nm*12], res.x[T*nm*12+1:end]

X = MPC.unbatch_states(x, 12)
p1 = plot(trans_states(flip_z(X))', labels=reshape(MPC.STATE_LABELS[MPC.TRANSLATIONAL_STATES], 1,6))
p2 = plot(reshape(u, 6, length(u) ÷ 6)', labels=reshape(["u$i" for i ∈ 1:6],1,6))
plot(p1, p2)







##

# u = repeat(basis(6,1)*10, nm*(T-1))
u = repeat(ones(6)*10, nm*(T-1))
_x = sys.A*repeat(x0, nm) + sys.B*u# - sys.Δ_nom




histogram(abs.(x - _x))

X1, X2 = MPC.unbatch_and_disjoint(_x, nm, T, 12)
