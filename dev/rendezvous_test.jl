using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Plots
using Clarabel

default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

begin
    ground = 10
    side = 6

    ylower = 0
    yupper = 4

    # ground = 15
    # side = 8

    γside = 0.5e-1
    γground = 0.5e-1

    constraints = MPC.RendezvousBarrier(h=ground, w=side, l_lower=ylower, l_upper=yupper,
                                          γg = γground, γs = γside)
    #constraints = [MPC.LinearConstraint(-basis(12, 3)*1, ground, γground)]

    #constraints = MPC.ElevatorShaft(h=6, γ=1.0)#γ=1e-1)
    modes = 0:1
    T = 30
    Δt = 10#0.05

    #u_bounds = (-Inf,Inf)
    #u_bounds = (0.1, 20.0)
    u_bounds = (-0.1, 0.1)

    num_modes = length(modes)

    nvec = [0.061, 0.090]
    models = [MPC.RendezvousModel(nvec[i]) for i in 1:num_modes]
    sys = MPC.BatchDynamics(models; T, Δt, u_bounds)

    ns = sys.inner_statedim
    nm = sys.inner_controldim
    x0 = zeros(ns)
    x0[1:3] .= [100, 3.8, 0]
    x_ref = zeros(ns)
    x_ref[1:3] .= [0, 0.5, 0]
    #x_ref[1:3] .= [5, -5, 10]

    #Qcustom = I(12) * 1.0e-1

    Qcustom = I(ns) * 1.0e-1
    Qcustom[1,1] = 50
    Qcustom[2,2] = 50

    Qcustom[3,3] = 50

    # Qcustom[3,3] = 1e1
    f = BarrierJuMPFormulator(
        sys,
        Clarabel.Optimizer;
        x_ref,
        Q=Qcustom,
        R=I(nm)* 1.0e-2,#1e-6,
        constraints,
        tol_feas = 1e-8,
        tol_infeas_abs = 1e-8,
        tol_infeas_rel = 1e-8, 
        tol_gap_abs = 1.0e-08,
        tol_gap_rel = 1.0e-08,
        verbose = true,
        max_iter= 50_000
    )
    model = JuMPModel(f, x0)

    #plot(unit_hist, lw=2)
    simtime = 40
    # failtime = floor(simtime/8)
    failtime = 2

    failmode = 2
    delaytime = 1
end

begin
    unit_planner = MPC.FTMPCPlanner(model, f, 1)
    #a, info = MPC.action_info(unit_planner, rand(12)*1e-2)

    unit_sim = Simulator(unit_planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
    unit_hist = simulate(unit_sim)
    #plot(unit_hist, lw=2)
    p = plot(unit_hist, Δt, side, ground)
end

begin
    consensus_planner = MPC.ConsensusSearchPlanner(model, f)
    consensus_sim = Simulator(consensus_planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
    function simcon(consensus_sim)
        return simulate(consensus_sim)
    end
    #consensus_hist = simulate(consensus_sim)
    consensus_hist = simcon(consensus_sim)
    plot(consensus_hist, lw=2, true)
    p = plot(consensus_hist, Δt, side, ground)

    #savefig(p, "consensus_hist.png")
    hists = [unit_hist, consensus_hist]
    totalplt = plot(hists, Δt, [side, ylower, yupper, ground], 
                    x_ref, Int(failtime), Int(delaytime)) |> display
    #plot(consensus_hist.consensus) |> display
end
begin
    animhist = unit_hist
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
            title="Max Consensus"
        )
    scatter3d!(plt, [x_ref[1]], [x_ref[2]], [-x_ref[3]], color="green")
    anim = @animate for i=1:maxanim
        scatter3d!(plt, [x[i]], [y[i]], [z[i]],
            xlims=(-7,7),
            ylims=(-7,7),
            zlims=(-15,5),
            camera = (40, 30),
            #color = i>partialtime ? "green" : "red"
            color = i>failtime ? "blue" : "red"
            )
        scatter3d!(plt, [x[i]], [y[i]], [zlims(plt)[1]], color="gray")
    end
    display(gif(anim, "anim_fps15.gif", fps=10))
end