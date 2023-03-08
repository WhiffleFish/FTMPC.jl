using ControlSystems
using Plots
using LinearAlgebra
using BarrierFTMPC
const MPC = BarrierFTMPC

model = LinearHexModel(0)
sys = model.ss
Δt = 0.1
dsys =  c2d(model.ss, Δt)
x0 = zeros(12)
x0[7] = 1.0
Ā, B̄ = batch_dynamics(dsys, 20)

x = Ā*x0

X = unbatch_states(x, 12)

x[[(i-1)*12 + 1 for i in 1:20]]

n_particles = 10
Ps = [zeros(12) for i ∈ 1:n_particles]
X0, A, B = joint_dynamics(Ps, sys.A, sys.B)
B
C = MPC.default_c(B)
D = MPC.default_d(B)

sys_big = ss(A,B,C,D)
x_ref = zeros(12)
x_ref[3] = -5
X_ref = reduce(vcat, Iterators.repeated(x_ref, n_particles))
Q = diagm(ones(12*n_particles))
R = diagm(ones(6))

L = lqr(sys_big, Q, R)
