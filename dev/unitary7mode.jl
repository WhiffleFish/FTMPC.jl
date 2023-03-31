using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using OSQP
using COSMO
using LinearAlgebra
using Plots

constraints = [
    #LinearConstraint(basis(12, 3)*1, 15, 1e-0),
    LinearConstraint(-basis(12, 3)*1, 3, 1e-1),
    LinearConstraint(basis(12, 2)*1, 1, 1e-0),
    LinearConstraint(-basis(12, 2)*1, 1, 1e-0),
    LinearConstraint(basis(12, 1)*1, 1, 1e-0),
    LinearConstraint(-basis(12, 1)*1, 1, 1e-0)
]

num_modes = 7
failures = 0:num_modes-1
# failures = [0,1]
T = 10
Δt = 0.05
u_bounds = (-Inf,Inf)
u_bounds = (.0,15.)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1] = 0
x_ref[2] = 0
x_ref[3] = 2

ws = [1,0,0,0,0,0,0]
#ws = [1,0]


Q_i = Matrix{Float64}(I(12))
Q_i[1,1] = 10.
#Q_i[3,3] = 10.
Q = [Q_i for i ∈ 1:num_modes] .* ws
R = [I(6)*0.01 for i ∈ 1:num_modes] .* ws

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    # P = (I(12)*1e-1,1),
    # Q = (I(6)*1e-6,1),
    P = Q,
    Q = R,
    constraints,
    # eps_prim_inf = 1e-3,
    # eps_abs = 1e-4,
    # eps_rel = 1e-4,
    verbose = false,
    # check_termination = 100,
    max_iter = 100_000
)

model = JuMPModel(f, x0)
#MPC.set_consensus_horizon(model, f, 4)
#= optimize!(model)

res = MPC.HexOSQPResults(f, model)
plot(pos_states(res.X[1])')
plot(pos_states(res.X[2])')
plot(pos_states(res.X[3])')

plot(res.U[1]')
plot(res.U[2]')
plot(res.U[3]') =#
simT = 40
imm = MPC.HexIMM(Δt=Δt)
planner = MPC.UnitaryConsensusPlanner(model, f)
sim = Simulator(imm, planner, x0=x0, T=simT)
hist,partialtime = simulate(sim)

display(plot(hist))
display(plot([hist.x'[:,[1,2]] -hist.x'[:,3]], label = ["x" "y" "z"]))
#plot(hist)
begin
    maxanim = size(hist.x)[2]
    xyz = [hist.x'[:,[1,2]] -hist.x'[:,3]]
    #plot([hist.x'[:,[1,2]] -hist.x'[:,3]], label = ["x" "y" "z"])
    x=xyz[:,1];y=xyz[:,2];z=xyz[:,3]
    zlow = extrema(z)[1]
    plt = scatter3d(
            xlims = extrema(x),
            ylims = extrema(y),
            zlims = extrema(z),
            xlabel = "X",
            ylabel = "Y",
            zlabel = "Z",
            #margin = 5mm,
            size = (800, 600),
            grid=true,
            legend=false,
            title="Robust-Maxh"
        )
    scatter3d!(plt, [x_ref[1]], [x_ref[2]], [-x_ref[3]], color="green")
    anim = @animate for i=1:maxanim
        scatter3d!(plt, [x[i]], [y[i]], [z[i]],
            xlims=(-2,10),
            ylims=(-1,1),
            zlims=(-10,11),
            camera = (40, 30),
            #color = i>partialtime ? "green" : "red"
            color = i>floor(simT/4) ? "blue" : "red"
            )
        scatter3d!(plt, [x[i]], [y[i]], [zlims(plt)[1]], color="gray")
    end
    display(gif(anim, "anim_fps15.gif", fps=10))
end


#plot(hist.x'[:,1:3])

#plot(hist.w[3,:])
#plot(hist.mode)

#plot(hist.w[2,:])

#= imm.T*basis(7,1)

model = MPC.set_objective_weights(planner, basis(7,4))
#MPC.set_consensus_horizon(planner, 10)
optimize!(model)

plot(hist)
hist.mode
hist.w[:,20]
MPC.weighted_sample(imm.T[:,1])

sys = imm.modes[1]
MPC.dstep(sys, zeros(12), ones(6)*10)
imm.u_noms[2]

imm.weights
 =#