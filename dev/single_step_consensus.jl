using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Plots
default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

constraints = MPC.ElevatorShaft(h=2, γ=1e-1)
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
    R=I(6)*1e-3,
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-4,
    eps_rel = 1e-4,
    verbose = false,
    max_iter= 50_000
)
model = JuMPModel(f, x0)

unit_planner = MPC.FTMPCPlanner(model, f, 1)
unit_sim = Simulator(MPC.HexIMM(;Δt), unit_planner, x0=x0, T=100, failure=MPC.FixedFailure(10,3;instant_update=true))
unit_hist = simulate(unit_sim)
plot(unit_hist)

consensus_planner = MPC.ConsensusSearchPlanner(model, f)
consensus_sim = Simulator(MPC.HexIMM(;Δt), consensus_planner, x0=x0, T=100, failure=MPC.FixedFailure(10,3;instant_update=true))
consensus_hist = simulate(consensus_sim)
plot(consensus_hist)

##
plot(consensus_hist.t, pos_states(consensus_hist.x)')
plot(consensus_hist.t, pos_states(consensus_hist.x)')
plot!(unit_hist.t[1:size(unit_hist.x,2)],pos_states(unit_hist.x)')

plot(consensus_hist.x[9,:])
plot!(unit_hist.x[9,:])

##
x = [-4.203591731230117e-12, 1.5069486592304665e-11, 0.4394728987947938, 2.881828239163578e-11, 9.752623246537374e-12, -2.1959617890069273e-12, -3.1180112935247936e-11, 1.0703823168356009e-10, 2.9251350911250324, -1.1432643852532679e-10, -2.6012220550684938e-11, -1.0836391228573638e-11]

MPC.set_initialstate(unit_planner, x)
optimize!(unit_planner.model)

MPC.action_info(unit_planner, x)


MPC.action(consensus_planner, last_x)





## extra step-wise MPC solution inspection
for i ∈ eachindex(consensus_hist.info)
    sleep(0.05)
    p = plot(consensus_hist.info[i][1])
    ylims!(p[1], (-2.2,0.1))
    ylims!(p[2], (-2,10))
    display(p)
end


for i ∈ eachindex(unit_hist.info)
    sleep(0.05)
    p = plot(unit_hist.info[i][1])
    ylims!(p[1], (-1.2,0.1))
    ylims!(p[2], (-2,10))
    display(p)
end

#= PLOT:
- 2 plots - pre-fail mode
    - Show MPC position for pre-fail mode pre-failure
    - Show MPC position for pre-fail mode post-failure
- 2 plots - post-fail mode
    - Show MPC position for post-fail mode pre-failure
    - Show MPC position for post-fail mode post-failure
=#

prefail_idx = 5
prefail_mode = 1
postfail_idx = 15
postfail_mode = 3

p1 = plot(
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][prefail_mode].X)))'
)
p2 = plot(
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][prefail_mode].X)))'
)
p3 = plot(
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][postfail_mode].X)))'
)
p4 = plot(
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][postfail_mode].X)))'
)

plot(p1,p2,p3,p4)

##


prefail_idx = 5
prefail_mode = 1
postfail_idx = 25
postfail_mode = 3

kwargs = (
    ylims = (-2.0,0.1),
    lw = 2
)
ts = 0:Δt:Δt*(T-1)

p1 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][prefail_mode].X)))';
    kwargs..., xticks=false
)
p2 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][prefail_mode].X)))';
    kwargs..., xticks=false, yticks=false
)
p3 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][postfail_mode].X)))';
    kwargs...
)
p4 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][postfail_mode].X)))';
    kwargs..., yticks = false
)

plot(p1,p2,p3,p4)
