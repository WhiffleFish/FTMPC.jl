using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using OSQP
using LinearAlgebra
using Plots

constraints = [
    LinearConstraint(basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 3)*0.5, 1, 1e-1)
]

failures = [0,1]
T = 20
Δt = 0.1
# u_bounds = (0.,10.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x0[9] = 10
x_ref = zeros(12)
x_ref[1:3] .= -5

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    P = I(12),#(I(12),1),
    Q = I(6)*1e-1,# (I(6)*1e-2,1),
    constraints
)

model = JuMPModel(f, x0)
res = MPC.HexOSQPResults(f, model)

p1 = plot(res.t,
    trans_states(flip_z(res.X[1]))',
    labels=permutedims(MPC.STATE_LABELS[MPC.TRANSLATIONAL_STATES]),
    lw=2
)
p2 = plot(res.t[1:end-1],res.U[1]', lw=2, labels=permutedims(["u$i" for i ∈ 1:6]))
plot(p1,p2)

##
using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using OSQP
using LinearAlgebra
using Plots

constraints = nothing

failures = [0,1]
T = 50
Δt = 0.1
# u_bounds = (0.,10.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    P = I(12),#(I(12),1),
    Q = I(6)*1e-2,# (I(6)*1e-2,1),
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-4,
    eps_rel = 1e-4,
    polish = true
)

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
