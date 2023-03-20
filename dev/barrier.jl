using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using OSQP
using LinearAlgebra
using Plots

A_constraint = zeros(12)
A_constraint[2] = -1
b_constraint = 1

failures = [0,1]
T = 100
Δt = 0.1
γ_constraint = 1e-1
# u_bounds = (0.,10.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[2] = 10
# x_ref[1:3] .= -5

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    P = I(12),#(I(12),1),
    Q = I(6)*1e-2,# (I(6)*1e-2,1),
    A_constraint,
    b_constraint,
    γ_constraint,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-4,
    eps_rel = 1e-4,
    polish = true
)

model = JuMPModel(f, x0)
MPC.set_initialstate(model, sys, zeros(12))
optimize!(model)
res = MPC.HexOSQPResults(f, model)

p1 = plot(res.t,
    trans_states(flip_z(res.X[1]))',
    labels=permutedims(MPC.STATE_LABELS[MPC.TRANSLATIONAL_STATES]),
    lw=2
)
p2 = plot(res.t[1:end-1],res.U[1]', lw=2, labels=permutedims(["u$i" for i ∈ 1:6]))
plot(p1,p2)

MPC.set_initialstate(model, sys, ones(12))

x = model[:x]
@objective(model, Min, 0.5*dot(x, f.P, x) + dot(f.q,x))

value.(model[:x])


##
barrier = f.barrier
x = value.(model[:x])

value.(model[:BARRIER])

value.(model[:BARRIER]) .≤ 0.5

@edit value(model[:BARRIER][2])

@edit JuMP._constraint_primal(model[:BARRIER][2],1)

max_barrier_violation = maximum(barrier.A*x .- barrier.ub)

##
using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using MadNLP
using LinearAlgebra
using Plots

A_constraint = zeros(12)
A_constraint[2] = -1
b_constraint = 1
γ_constraint = 1e-1

failures = [0]
T = 100
Δt = 0.1

# u_bounds = (0.,10.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[2] = 10
# x_ref[1:3] .= -5

f = BarrierJuMPFormulator(
    sys,
    MadNLP.Optimizer;
    x_ref, Q=I(6)*1e-2,
    A_constraint,
    b_constraint,
    γ_constraint
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


get_kwargs(;kwargs...) = kwargs
kwargs = get_kwargs(a=1, b=2)
convert_kwargs(kwargs::Base.Pairs) = Tuple(string(a)=>b for (a,b) ∈ kwargs)

convert_kwargs(kwargs)


##
using ControlSystems
s = LinearHexModel()
ss = c2d(s.ss, 0.1)
step(ss, rand(6))
