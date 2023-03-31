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

ground = 6
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
failmode = 1
failures = [0,failmode]
#failures = 0:num_modes-1
# failures = [0,1]
T = 10
Δt = 0.05
u_bounds = (-Inf,Inf)
u_bounds = (.0,15.)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1] = -0.4
x_ref[2] = 0.4
x_ref[3] = 5

#ws = [1,0,0,0,0,0,0]
ws = [1,0]


Q_i = Matrix{Float64}(I(12))
Q_i[1,1] = Q_i[2,2] = 20.
#Q_i[3,3] = 10.
Q = [Q_i for i ∈ 1:num_modes] .* ws
R = [I(6)*0.01 for i ∈ 1:num_modes] .* ws

f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    P = Q,
    Q = R,
    eps_abs = 1e-4,
    eps_rel = 1e-4,
    constraints,

    verbose = false,
    max_iter = 100_000
)

model = JuMPModel(f, x0)

simT = 80
failtime = floor(simT/4)
delaytime = 3
imm = MPC.HexIMM(Δt=Δt)

planner = MPC.UnitaryConsensusPlanner(model, f)
sim = Simulator(imm, planner, x0=x0, T=simT)
unithist,partialtime = simulate(sim, failmode, failtime, delaytime)


nominalplanner = MPC.ConsensusSearchPlanner(model, f)
nominalsim = Simulator(imm, nominalplanner, x0=x0, T=simT)
nominalhist,partialtime = simulate(nominalsim, failmode, failtime, delaytime)



num_modes = 1
failures = [0]
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
ws = [1]
Q = [Q_i for i ∈ 1:num_modes] .* ws
R = [I(6)*0.01 for i ∈ 1:num_modes] .* ws
fnr = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    P = Q,
    Q = R,
    constraints,
    eps_abs = 1e-4,
    eps_rel = 1e-4,
    verbose = false,
    max_iter = 100_000
)

modelnr = JuMPModel(fnr, x0)
nonrobustplanner = MPC.NonRobustPlanner(modelnr, fnr)
nonrobustsim = Simulator(imm, nonrobustplanner, x0=x0, T=simT)
nonrobusthist = simulate(nonrobustsim, failmode, failtime, delaytime)

print("Done")


# PLOTTING 
firstplt = plot(unithist, Δt, side, ground)
#secondplt = plot(nominalhist, Δt, side, ground)
nonrobplt = plot(nonrobusthist, Δt, side, ground)

crashind = findfirst(x -> x>1, nonrobusthist.x[1,:])
hists = [nonrobusthist, unithist, nominalhist]
totalplt = plot(hists, Δt, side, ground, x_ref, Int(failtime), Int(delaytime))
savefig(totalplt, joinpath(@__DIR__, "../figs/comparelowtol.svg"))

display(totalplt)
#display(firstplt)
#display(secondplt)

# Save to matfile
#= begin
    states = getfield.(hists, :x)
    file = matopen("matfile.mat", "w")
    write(file, "states", states)
    close(file)
end =#

# Animation
begin
    animhist = nominalhist
    xyz = [animhist.x'[:,[1,2]] -animhist.x'[:,3]]
    x=xyz[:,1];y=xyz[:,2];z=xyz[:,3]
    tsim = length(x)
    tvec = Δt:Δt:tsim*Δt
    maxanim = size(animhist.x)[2]
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
            color = i>failtime ? "blue" : "red"
            )
        scatter3d!(plt, [x[i]], [y[i]], [zlims(plt)[1]], color="gray")
    end
    display(gif(anim, "anim_fps15.gif", fps=10))
end

