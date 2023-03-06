using ControlSystems
using Plots
using LinearAlgebra
using BarrierFTMPC
const MPC = BarrierFTMPC


model = LinearHexModel()
x0 = zeros(12)
x_ref = zeros(12)
x_ref[3] = -5

Q_v = zeros(12)
Q_v[3] = 1.0
Q = diagm(Q_v)
R = zeros(6)

Q = I
R = I
L = lqr(model.ss, Q, R)

u(x,t) = -L*(x - x_ref)

t = 0:0.1:10
lsim(model.ss, u, t, x0=x0)
y, t, x, uout = lsim(model.ss,u,t,x0=x0)

plot(t, trans_states(flip_z(x))', labels=MPC.STATE_LABELS)
plot(t, uout')

plot(x[7:end,:]')
