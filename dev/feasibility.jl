begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    using BarrierFTMPC
    const MPC = BarrierFTMPC
    using JuMP
    using LinearAlgebra
    Pkg.activate(@__DIR__)
    using OSQP
    using COSMO
end
using Plots

#= LP Solvers
- Clp
- CDCS
- CDD
- Clarabel
- COSMO
=#

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
    P = (I(12),1),#(I(12),1),
    Q = (I(6)*1e-2,1),# (I(6)*1e-2,1),
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-5,
    eps_rel = 1e-5,
    verbose = true
)

model = JuMPModel(f, x0)

f_model = copy(model)
@objective(f_model, Min, 0)
set_optimizer(f_model, OSQP.Optimizer)
optimize!(f_model)
set_optimizer(f_model, COSMO.Optimizer)
optimize!(f_model)


optimize!(f_model)

set_objective_function(f_model, nothing)


planner = MPC.ConsensusSearchPlanner(model, f)
sim = Simulator(LinearHexModel(0), planner, x0=x0, T=50)
hist = simulate(sim)

plot(hist)
