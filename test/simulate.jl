failures = [0,1]
T = 50
Δt = 0.1
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5
ss = LinearHexModel(0)

f = JuMPFormulator(sys, OSQP.Optimizer;x_ref,Q=I(6)*1e-2, verbose=false)
model = JuMPModel(f, x0)
optimize!(model) # warm start

planner = FTMPCPlanner(model, f)
sim = MPC.Simulator(LinearHexModel(0), planner, x0=zeros(12), T=50)
hist = simulate(sim)

@test hist.x[:, end] ≈ x_ref atol = 0.25

## with barriers

constraints = [
    LinearConstraint(basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 3)*0.5, 1, 1e-1)
]

failures = [0,1]
T = 30
Δt = 0.1

u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = basis(12,9)*3
x_ref = zeros(12)
x_ref[1:3] .= -10

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    Q=I(6)*1e-1,
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-5,
    eps_rel = 1e-5,
    verbose = false,
    max_iter= 5_000
    # time_limit = Δt*1.5
)
model = JuMPModel(f, x0)
optimize!(model)

res = MPC.HexOSQPResults(f, model)

planner = FTMPCPlanner(model, f)
sim = MPC.Simulator(LinearHexModel(0), planner, x0=x0, T=100)
hist = simulate(sim)

@test hist.x[1:2,end] ≈ x_ref[1:2] atol = 1e-1
@test hist.x[4:end,end] ≈ x_ref[4:end] atol = 1e-1
@test hist.x[3,end] ≈ -1.0 atol = 1e-1 # close to barrier

@test MPC.max_barrier_violation(sim, hist) ≤ 0.1
