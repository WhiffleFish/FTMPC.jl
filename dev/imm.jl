using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using OSQP
using COSMO
using LinearAlgebra
using Plots

constraints = [
    LinearConstraint(basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(basis(12, 2)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 2)*1, 1, 1e-1),
    # LinearConstraint(basis(12, 1)*1, 1, 1e-1),
    # LinearConstraint(-basis(12, 1)*1, 1, 1e-1)
]

failures = 0:6
T = 20
Δt = 0.1
u_bounds = (-Inf,Inf)
u_bounds = (.1,15.)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= 5

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    P = (I(12)*1e-3,1),
    Q = (I(6)*1e-5,1),
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-3,
    eps_rel = 1e-3,
    verbose = true,
    check_termination = 100,
    rho = 1000.,
    max_iter = 100_000
)

model = JuMPModel(f, x0)
planner = MPC.ConsensusSearchPlanner(model, f)
sim = Simulator(MPC.HexIMM(), planner, x0=x0, T=30)
hist = simulate(sim)

plot(hist)
plot(hist.x'[:,1:3])
plot(hist.u')

plot(hist.w[3,:])
plot(hist.mode)

plot(hist.w[2,:])

imm.T*basis(7,1)

model = MPC.set_objective_weights(planner, basis(7,4))
MPC.set_consensus_horizon(planner, 10)
optimize!(model)

plot(hist)
hist.mode
hist.w[:,20]
MPC.weighted_sample(imm.T[:,1])

sys = imm.modes[1]
MPC.dstep(sys, zeros(12), ones(6)*10)
imm.u_noms[2]

imm.weights
