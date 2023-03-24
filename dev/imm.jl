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
    # P = I(12)*1e-4,
    P = (I(12),1),
    # Q = I(6)*1e-5,
    Q = (I(6)*1e-2,1),
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-4,
    eps_rel = 1e-4,
    verbose = true,
    max_iter = 100_000
)

model = JuMPModel(f, x0)
planner = MPC.ConsensusSearchPlanner(model, f)
imm = MPC.HexIMM()
sim = Simulator(imm, planner, x0=x0, T=50)
hist = simulate(sim)

plot(hist)

plot(hist.x[3,:])

hist.w
hist.mode



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


##

(; weights, modes, u_noms, T, obs_dist) = imm
true_mode = 2
x = zeros(12)
u = ones(6)
δu_true = u - u_noms[true_mode]
xp = MPC.dstep(modes[true_mode], x, δu_true)
y = rand(obs_dist(xp))
pdf(obs_dist(xp), y)

weights .= basis(7,1)
weights = T*weights # predictor
_mode = 1
sys = modes[_mode]
u_nom = u_noms[_mode]
δu = u - u_nom
xp_guess = MPC.dstep(sys, x, δu)
@show pdf(obs_dist(xp_guess), y)
weights[mode] *= pdf(obs_dist(xp), y) # corrector

for mode in eachindex(modes)
    sys = modes[mode]
    u_nom = u_noms[mode]
    δu = u - u_nom
    xp = dstep(sys, x, δu)
    @show xp
    @show pdf(obs_dist(xp), y)
    weights[mode] *= pdf(obs_dist(xp), y) # corrector
end
