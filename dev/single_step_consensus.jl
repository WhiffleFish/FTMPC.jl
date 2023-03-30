using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using COSMO
using LinearAlgebra
using Plots
default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

constraints = MPC.ElevatorShaft(h=1, γ = 1e-1)

failures = 0:6
T = 10
Δt = 0.05

u_bounds = (-15.,15.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = basis(12,3)*(2)


f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    Q=I(12)*1e-1,
    R=I(6)*1e-5,
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-4,
    eps_rel = 1e-4,
    verbose = false,
    max_iter= 50_000
)
model = JuMPModel(f, x0)
planner = MPC.FTMPCPlanner(model, f, 1)
sim = Simulator(MPC.HexIMM(;Δt), planner, x0=x0, T=100, failure=MPC.FixedFailure(10,3;instant_update=true))
hist = simulate(sim)
plot(hist)

planner = MPC.ConsensusSearchPlanner(model, f)
sim = Simulator(MPC.HexIMM(;Δt), planner, x0=x0, T=100, failure=MPC.FixedFailure(10,3;instant_update=true))
hist = simulate(sim)
plot(hist)


## extra step-wise MPC solution inspection
for i ∈ eachindex(hist.info)
    sleep(0.05)
    p = plot(hist.info[i][1])
    ylims!(p[1], (-1.2,0.1))
    ylims!(p[2], (-2,10))
    display(p)
end

