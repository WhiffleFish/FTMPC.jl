using ControlSystems
using Plots
using LinearAlgebra
using BarrierFTMPC
const MPC = BarrierFTMPC


begin
    model = LinearHexModel(2)
    x0 = zeros(12)
    x_ref = zeros(12)
    x_ref[3] = -5
    Q = I; R = I
    L = lqr(model.ss, Q, R)
    u(x,t) = -L*(x - x_ref)
    t = 0:0.1:10
    lsim(model.ss, u, t, x0=x0)
    y, t, x, uout = lsim(model.ss,u,t,x0=x0)
    p1 = plot(t, trans_states(flip_z(x))', labels=MPC.STATE_LABELS)
    p2 = plot(t, (uout .+ model.u)', labels=reshape(["u$i" for i ∈ 1:6],1,6), lw=2)
    plot(p1,p2)
end

plot(t, uout', labels=reshape(["δu$i" for i ∈ 1:6],1,6), lw=2)
plot(x[7:end,:]')

(uout .+ model.u)[:,end]
uout[:,end]
