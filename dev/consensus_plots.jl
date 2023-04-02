using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Plots
default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

constraints = MPC.ElevatorShaft(h=6, γ=1e-1)
failures = 0:6
T = 10
Δt = 0.05

u_bounds = (-15.,15.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= [-0.4, 0.4, 5]

q_vec = zeros(12)
q_vec[1:3] .= 1.
Q = diagm(q_vec)


f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    Q=I(12),
    R=I(6)*1e-3,
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-5,
    eps_rel = 1e-5,
    verbose = false,
    max_iter= 50_000
)
model = JuMPModel(f, x0)

unit_planner = MPC.FTMPCPlanner(model, f, 1)
unit_sim = Simulator(MPC.HexIMM(;Δt), unit_planner, x0=x0, T=100, failure=MPC.FixedFailure(10,3;instant_update=true))
unit_hist = simulate(unit_sim)
plot(unit_hist, lw=2)

consensus_planner = MPC.ConsensusSearchPlanner(model, f)
consensus_sim = Simulator(MPC.HexIMM(;Δt), consensus_planner, x0=x0, T=100, failure=MPC.FixedFailure(20,3;instant_update=true))
consensus_hist = simulate(consensus_sim)
plot(consensus_hist)

using LaTeXStrings
p = plot(
    consensus_hist.t, consensus_hist.consensus, 
    ylabel=L"h^*", xlabel="Time (s)", 
    lw=2, title="Consensus Horizon History",
    xlims=(0,5),
    grid=true,
    yticks = 0:10
)
savefig(p, "consensus_horizon_history.pdf")

using FileIO
using JLD2
save("hist.jld2", Dict("hist"=>consensus_hist))

##
plot(consensus_hist.t, pos_states(consensus_hist.x)')
plot!(unit_hist.t[1:size(unit_hist.x,2)],pos_states(unit_hist.x)')

## extra step-wise MPC solution inspection
for i ∈ eachindex(consensus_hist.info)
    sleep(0.05)
    p = plot(consensus_hist.info[i][3])
    ylims!(p[1], (-6.2,0.1))
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

Ts = 0:Δt:Δt*(consensus_sim.T-1)

prefail_idx = 5
prefail_time = Ts[prefail_idx]
prefail_mode = 1
postfail_idx = 51
postfail_time = Ts[postfail_idx]
postfail_mode = 3

kwargs = (
    ylims = (-6.0,0.5),
    lw = 2,
    # grid = true
)
ts = 0:Δt:Δt*(T-1)

p1 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][prefail_mode].X)))';
    kwargs..., xticks=false, ylabel = "Nominal Mode", title = "Pre-failure",
    labels = ["x" "y" "z"]
)
p2 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][prefail_mode].X)))';
    kwargs..., xticks=false, yticks=false, title = "Post-failure"
)
p3 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][postfail_mode].X)))';
    kwargs..., ylabel = "Failing Mode"
)
p4 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][postfail_mode].X)))';
    kwargs..., yticks = false
)

plot(p1,p2,p3,p4)

##

prefail_idx = 5
prefail_time = Ts[prefail_idx]
prefail_mode = 1
fail_idx = 21
fail_time = Ts[fail_idx]
postfail_idx = 51
postfail_time = Ts[postfail_idx]
postfail_mode = 3

kwargs = (
    ylims = (-6.0,0.5),
    lw = 2,
    grid = true
)
ts = 0:Δt:Δt*(T-1)

p1 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][prefail_mode].X)))';
    kwargs..., x_formatter = _->"", ylabel = "Nominal Mode", title = "Pre-failure",
    labels = ["x" "y" "z"]
)
p2 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[fail_idx][prefail_mode].X)))';
    kwargs..., x_formatter = _->"", y_formatter = _->"", title = "Failure"
)
p3 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][prefail_mode].X)))';
    kwargs..., x_formatter = _->"", y_formatter = _->"", title = "Post-failure"
)
p4 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][postfail_mode].X)))';
    kwargs..., ylabel = "Failing Mode"
)
p5 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[fail_idx][postfail_mode].X)))';
    kwargs..., y_formatter = _->""
)
p6 = plot(ts,
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][postfail_mode].X)))';
    kwargs..., y_formatter = _->""
)
_plots = (p1,p2,p3,p4,p5,p6)
for (i,p) ∈ enumerate(_plots)
    hline!(p, [-x_ref[3]], c=:firebrick, ls=:dash, label=isone(i) ? "ref" : "", lw=2)
end

plot(_plots...)

##




kwargs = (
    ylims = (-6.0,0.5),
    lw = 2,
    grid = true,
    titlefontsize = 10,
    xguidefontsize = 10,
    yguidefontsize = 10,
)

p1 = plot(
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][prefail_mode].X)))';
    kwargs..., x_formatter = _->"", ylabel = "Nominal Mode", title = "Pre-failure (t=$prefail_time s)",
    labels = ["x" "y" "z"]
)
p2 = plot(
    pos_states(flip_z(only(consensus_hist.info[fail_idx][prefail_mode].X)))';
    kwargs..., x_formatter = _->"", y_formatter = _->"", title = "Failure (t=$fail_time s)"
)
p3 = plot(
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][prefail_mode].X)))';
    kwargs..., x_formatter = _->"", y_formatter = _->"", title = "Post-failure  (t=$postfail_time s)"
)
p4 = plot(
    pos_states(flip_z(only(consensus_hist.info[prefail_idx][postfail_mode].X)))';
    kwargs..., ylabel = "Failing Mode", xlabel = "Planning Horizon"
)
p5 = plot(
    pos_states(flip_z(only(consensus_hist.info[fail_idx][postfail_mode].X)))';
    kwargs..., y_formatter = _->"", xlabel = "Planning Horizon"
)
p6 = plot(
    pos_states(flip_z(only(consensus_hist.info[postfail_idx][postfail_mode].X)))';
    kwargs..., y_formatter = _->"", xlabel = "Planning Horizon"
)
_plots = (p1,p2,p3,p4,p5,p6)
for (i,p) ∈ enumerate(_plots)
    hline!(p, [-x_ref[3]], c=:firebrick, ls=:dash, label=isone(i) ? "ref" : "", lw=2)
end

p = plot(_plots...)
savefig(p, "optimizer_trajectories.pdf")

##
kwargs = (
    ylims = (-6.0,0.5),
    lw = 2,
    grid = true,
    titlefontsize = 10,
    xguidefontsize = 10,
    yguidefontsize = 10,
)

Ts = 0:Δt:Δt*(consensus_sim.T-1)
prefail_idx = 5
prefail_time = Ts[prefail_idx]
prefail_mode = 1
fail_idx = 21
fail_time = Ts[fail_idx]
postfail_idx = 51
postfail_time = Ts[postfail_idx]
postfail_mode = 3

process_info(consensus_hist.info[postfail_idx][prefail_mode])
process_info(info) = flip_z(only(info.X))[3,:]
p1 = plot(
    process_info(consensus_hist.info[prefail_idx][prefail_mode]);
    kwargs..., x_formatter = _->"", ylabel = "Nominal Mode", title = "Pre-failure (t=$prefail_time s)",
    label = "z"
)
p2 = plot(
    process_info(consensus_hist.info[fail_idx][prefail_mode]);
    kwargs..., x_formatter = _->"", y_formatter = _->"", title = "Failure (t=$fail_time s)"
)
p3 = plot(
    process_info(consensus_hist.info[postfail_idx][prefail_mode]);
    kwargs..., x_formatter = _->"", y_formatter = _->"", title = "Post-failure  (t=$postfail_time s)"
)
p4 = plot(
    process_info(consensus_hist.info[prefail_idx][postfail_mode]);
    kwargs..., ylabel = "Mode 2", xlabel = "Planning Horizon"
)
p5 = plot(
    process_info(consensus_hist.info[fail_idx][postfail_mode]);
    kwargs..., y_formatter = _->"", xlabel = "Planning Horizon"
)
p6 = plot(
    process_info(consensus_hist.info[postfail_idx][postfail_mode]);
    kwargs..., y_formatter = _->"", xlabel = "Planning Horizon"
)
_plots = (p1,p2,p3,p4,p5,p6)
for (i,p) ∈ enumerate(_plots)
    hline!(p, [-x_ref[3]], c=:firebrick, ls=:dash, label=isone(i) ? "ref" : "", lw=2)
end

p = plot(_plots...)
savefig(p, "optimizer_trajectories.pdf")
