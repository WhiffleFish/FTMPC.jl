constraints = [
    LinearConstraint(basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 3)*0.5, 1, 1e-1)
]

failures = [0,1]
T = 100
Δt = 0.1
# u_bounds = (0.,10.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[2] = 10

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
    max_iter = 10_000
)

model = JuMPModel(f, x0)
optimize!(model)
res = MPC.OSQPResults(f, model)

@test MPC.max_barrier_violation(f, model) ≤ 0.1

opt_infeasibility = maximum(value.(model[:BARRIER]) .- normalized_rhs.(model[:BARRIER]))
@test opt_infeasibility ≤ 1e-2
