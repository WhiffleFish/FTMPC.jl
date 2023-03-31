using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using OSQP
using COSMO
using LinearAlgebra
using Plots
default(grid=false, framestyle=:box, fontfamily="Computer Modern")
theme(:wong)

ground = 3
side = 1
γside = 1e-0
γground = 1e-1

constraints = MPC.ElevatorShaft(h=ground, γg = γground, γs = γside)

#= constraints = [
    #LinearConstraint(basis(12, 3)*1, 15, 1e-0),
    LinearConstraint(-basis(12, 3)*1, ground, γground),
    LinearConstraint(basis(12, 2)*1, side, γside),
    LinearConstraint(-basis(12, 2)*1, side, γside),
    LinearConstraint(basis(12, 1)*1, side, γside),
    LinearConstraint(-basis(12, 1)*1, side, γside)
]
 =#
num_modes = 2
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

#ws = [1,0,0,0,0,0,0]
ws = [1,0]


Q_i = Matrix{Float64}(I(12))
Q_i[1,1] = 10.
#Q_i[3,3] = 10.
Q = [Q_i for i ∈ 1:num_modes] .* ws
R = [I(6)*0.01 for i ∈ 1:num_modes] .* ws

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    P = Q,
    Q = R,

    constraints,

    verbose = false,
    max_iter = 100_000
)

model = JuMPModel(f, x0)

simT = 40
imm = MPC.HexIMM(Δt=Δt)
#= planner = MPC.UnitaryConsensusPlanner(model, f)
sim = Simulator(imm, planner, x0=x0, T=simT)
hist,partialtime = simulate(sim)

nominalplanner = MPC.ConsensusSearchPlanner(model, f)
nominalsim = Simulator(imm, nominalplanner, x0=x0, T=simT)
nominalhist,partialtime = simulate(nominalsim)
 =#
nonrobustplanner = MPC.NonRobustPlanner(model, f)
nonrobustsim = Simulator(imm, nonrobustplanner, x0=x0, T=simT)
nonrobusthist,partialtime = simulate(nonrobustsim)

print("Done")

# PLOTTING 
#firstplt = plot(hist, Δt, side, ground)
#secondplt = plot(nominalhist, Δt, side, ground)
nonrobplt = plot(nonrobusthist, Δt, side, ground)
#hists = [hist, nominalhist, nonrobusthist]
#totalplt = plot(hists, Δt, side, ground)
#display(totalplt)
#display(firstplt)
#display(secondplt)

begin
    xyz = [hist.x'[:,[1,2]] -hist.x'[:,3]]
    x=xyz[:,1];y=xyz[:,2];z=xyz[:,3]
    tsim = length(x)
    tvec = Δt:Δt:tsim*Δt
    maxanim = size(hist.x)[2]
    x=xyz[:,1];y=xyz[:,2];z=xyz[:,3]
    zlow = extrema(z)[1]
    plt = scatter3d(
            xlims = extrema(x),
            ylims = extrema(y),
            zlims = extrema(z),
            xlabel = "X",
            ylabel = "Y",
            zlabel = "Z",
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

