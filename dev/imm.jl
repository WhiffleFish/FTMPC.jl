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
    LinearConstraint(basis(12, 1)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 1)*1, 1, 1e-1)
]

failures = 0:6
# failures = [0,1]
T = 10
Δt = 0.01
# u_bounds = (-Inf,Inf)
u_bounds = (-15.,15.)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1] = 5
x_ref[2] = 10
x_ref[3] = -10

ws = [1,0,0,0,0,0,0]
ws = [1,0.1,0.1,0.1,0.1,0.1,0.1]

state_weights = zeros(12)
state_weights[1:3] .= 1

Q_i = diagm(state_weights)

# Q_i = Matrix{Float64}(I(12))
# Q_i[3,3] = 10.
# Q = [Q_i for i ∈ 1:7] .* ws
R = [I(6)*1e-10 for i ∈ 1:7] .* ws

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    # P = (I(12)*1e-1,1),
    # Q = (I(6)*1e-6,1),
    Q = Q_i,
    R = R,
    constraints,
    # eps_prim_inf = 1e-3,
    eps_abs = 1e-5,
    eps_rel = 1e-5,
    verbose = false,
    # check_termination = 100,
    max_iter = 100_000
)

model = JuMPModel(f, x0)
MPC.set_consensus_horizon(model, f, 1)
optimize!(model)

res = MPC.HexOSQPResults(f, model)
p1 = plot(pos_states(res.X[1])')
p2 = plot(pos_states(res.X[2])')
p3 = plot(pos_states(res.X[3])')

plot((plot(pos_states(res.X[i])') for i ∈ 1:6)...)

plot(res.U[1]')
plot(res.U[2]')
plot(res.U[3]')

planner = MPC.FTMPCPlanner(model, f, 2)
planner = MPC.ConsensusSearchPlanner(model, f)
sim = Simulator(MPC.HexIMM(;Δt), planner, x0=x0, T=100)
hist = simulate(sim)

plot(hist)
hist.x
plot(hist.x[1:3,:]')
plot(hist.u')

plot(trans_states(hist.info[10].X[6])')
# cost is somehow negative???

MPC.set_consensus_horizon(model, f, 7)
optimize!(model)

