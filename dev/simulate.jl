using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Plots
default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

failures = [0,1]
T = 50
Δt = 0.1
# u_bounds = (0.,12.)
# u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5
ss = LinearHexModel(0)

f = JuMPFormulator(sys, OSQP.Optimizer;x_ref,Q=I(6)*1e-2, polish=false, verbose=false)
model = JuMPModel(f, x0)
optimize!(model) # warm start

planner = FTMPCPlanner(model, f)
sim = MPC.Simulator(LinearHexModel(0), planner, x0=zeros(12), T=50)
hist = simulate(sim)

plot(hist, lw=2)

## with barriers

constraints = [
    LinearConstraint(basis(12, 3)*1, 1, 1e-1),
    LinearConstraint(-basis(12, 3)*0.5, 1, 1e-1)
]

failures = [0,1]
T = 30
Δt = 0.1

u_bounds = (-10.,10.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = basis(12,9)*3 # zeros(12)
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

planner = FTMPCPlanner(model, f)
sim = Simulator(LinearHexModel(0), planner, x0=x0, T=100)
hist = simulate(sim)

plot(hist, lw=2)


MPC.max_barrier_violation(sim,hist)
maximum(f.barrier.A*x .- f.barrier.ub)
f.barrier.A
ub = first(f.barrier.ub)
maximum(dot(f.barrier.A[1,1:12],x)-ub for x ∈ eachcol(hist.x))

f.barrier.A*value.(model[:x])


MPC.set_objective_weights(model, f, [1,0])
MPC.set_objective_weights(model, f, [1,0])

c = model[:CONSENSUS][10]


6*nm*(T-1)
