using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Plots
using Clarabel

default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")

begin
    ground = 6
    side = 1

    # ground = 15
    # side = 8

    γside = 0.5e-1
    γground = 0.5e-1

    constraints = MPC.ElevatorShaft(h=ground, w=side, l=side, γg = γground, γs = γside)
    #constraints = [MPC.LinearConstraint(-basis(12, 3)*1, ground, γground)]

    #constraints = MPC.ElevatorShaft(h=6, γ=1.0)#γ=1e-1)
    failures = 0:2
    T = 10
    Δt = 0.05#0.05

    #u_bounds = (-Inf,Inf)
    #u_bounds = (0.1, 20.0)
    u_bounds = (0.1, 20.0)

    nm = length(failures)

    models = [MPC.HexCTLinearModel(failure) for failure in failures]
    sys = MPC.BatchDynamics(models; T, Δt, u_bounds)

    x0 = zeros(12)
    x_ref = zeros(12)
    x_ref[1:3] .= [-0.7, 0.7, 5]
    #x_ref[1:3] .= [5, -5, 10]

    #Qcustom = I(12) * 1.0e-1

    Qcustom = I(12) * 1.0e-1
    Qcustom[1,1] = 50
    Qcustom[2,2] = 50

    Qcustom[3,3] = 50

    # Qcustom[3,3] = 1e1
    f = BarrierJuMPFormulator(
        sys,
        Clarabel.Optimizer;
        x_ref,
        Q=Qcustom,
        R=I(6)* 1.0e-2,#1e-6,
        constraints,
        tol_feas = 1e-8,
        tol_infeas_abs = 1e-8,
        tol_infeas_rel = 1e-8, 
        tol_gap_abs = 1.0e-08,
        tol_gap_rel = 1.0e-08,
        verbose = false,
        max_iter= 50_000
    )
    model = JuMPModel(f, x0)

    #plot(unit_hist, lw=2)
    simtime = 80
    failtime = floor(simtime/4)
    failmode = 2
    delaytime = 3
end

begin "Unit"
    unit_planner = MPC.FTMPCPlanner(model, f, 1)
    #a, info = MPC.action_info(unit_planner, rand(12)*1e-2)

    unit_sim = Simulator(unit_planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
    unit_hist = simulate(unit_sim)
    plot(unit_hist, lw=2)
    p = plot(unit_hist, Δt, side, ground)
end

begin "Full"
    full_planner = MPC.FTMPCPlanner(model, f, T-1)
    full_sim = Simulator(full_planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
    full_hist = simulate(full_sim)
    #plot(full_hist, lw=2)
    p = plot(full_hist, Δt, side, ground)
end

begin "Nonrobust"
    nonrobust_planner = MPC.FTMPCPlanner(model, f, 0)
    nonrobust_sim = Simulator(nonrobust_planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
    nonrobust_hist = simulate(nonrobust_sim)
    p = plot(nonrobust_hist, Δt, side, ground)
end

begin "Consensus"
    consensus_planner = MPC.ConsensusSearchPlanner(model, f)
    consensus_sim = Simulator(consensus_planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
    consensus_hist = simulate(consensus_sim)
    plot(consensus_hist, lw=2)
    p = plot(consensus_hist, Δt, side, ground)

    hists = [consensus_hist, unit_hist, nonrobust_hist, full_hist]
    x_ref[3] = -x_ref[3]
    totalplt = plot(hists, Δt, side, ground, x_ref, Int(failtime), Int(delaytime)) |> display
    #plot(consensus_hist.consensus) |> display
end

begin
    function draw_rectangle!(plt, vertices; color="blue", alpha=0.5)
        # Draw the four edges of the rectangle
        plot!(plt, [vertices[1][1], vertices[2][1]], [vertices[1][2], vertices[2][2]], [vertices[1][3], vertices[2][3]], color=color, alpha=alpha, legend=false) # Edge 1
        plot!(plt, [vertices[2][1], vertices[3][1]], [vertices[2][2], vertices[3][2]], [vertices[2][3], vertices[3][3]], color=color, alpha=alpha, legend=false) # Edge 2
        plot!(plt, [vertices[3][1], vertices[end][1]], [vertices[3][2], vertices[end][2]], [vertices[3][3], vertices[end][3]], color=color, alpha=alpha, legend=false) # Edge 3
        plot!(plt, [vertices[end][1], vertices[1][1]], [vertices[end][2], vertices[1][2]], [vertices[end][3], vertices[1][3]], color=color, alpha=alpha, legend=false) # Edge 4
    end
    
    hist_colors = [:cyan3, :maroon3, :blue, :gold2]
    rectangles_vertices = [
        [(-1, -1, -6), (1, -1, -6), (1, -1, 0), (-1, -1, 0)],
        [(1, -1, -6), (1, 1, -6), (1, 1, 0), (1, -1, 0)],
        [(-1, 1, -6), (1, 1, -6), (1, 1, 0), (-1, 1, 0)],
        [(-1, -1, -6), (-1, 1, -6), (-1, 1, 0), (-1, -1, 0)],
        [(-1, -1, -6), (1, -1, -6), (1, 1, -6), (-1, 1, -6)],
    ]

    animhist = consensus_hist
    xyz = [animhist.x'[:,[1,2]] -animhist.x'[:,3]]
    x_con=xyz[:,1];y_con=xyz[:,2];z_con=xyz[:,3]
    feasvec_con = [info.feas for info ∈ unit_hist.info]

    xyz_full = [full_hist.x'[:,[1,2]] -full_hist.x'[:,3]]
    x_full=xyz_full[:,1];y_full=xyz_full[:,2];z_full=xyz_full[:,3]

    xyz_unit = [unit_hist.x'[:,[1,2]] -unit_hist.x'[:,3]]
    x_unit=xyz_unit[:,1];y_unit=xyz_unit[:,2];z_unit=xyz_unit[:,3]

    xyz_nonrobust = [nonrobust_hist.x'[:,[1,2]] -nonrobust_hist.x'[:,3]]
    x_nonrobust=xyz_nonrobust[:,1];y_nonrobust=xyz_nonrobust[:,2];z_nonrobust=xyz_nonrobust[:,3]

    xvec = [x_con, x_unit, x_nonrobust, x_full]
    yvec = [y_con, y_unit, y_nonrobust, y_full]
    zvec = [z_con, z_unit, z_nonrobust, z_full]

    feasvec = [[info.feas for info ∈ hists[i].info] for i ∈ eachindex(hists)]

    tsim = length(x)
    tvec = Δt:Δt:tsim*Δt
    maxanim = size(animhist.x)[2]
    x=xyz[:,1];y=xyz[:,2];z=xyz[:,3]
    zlow = extrema(z)[1]

    pltvec = []
    titles = ["Maximize Consensus (Ours)", "First-Step Consensus", 
            "Non-Robust", "Full-Step Consensus"]
    for i ∈ eachindex(hists)
        push!(pltvec, scatter3d(
            xlims = extrema(x_con),
            ylims = extrema(y_con),
            zlims = [-6,2],
            xlabel = "X",
            ylabel = "Y",
            zlabel = "Z",
            grid=true,
            legend=false,
            titlefontsize=9,
            title=titles[i]
        ))
        for vertices in rectangles_vertices
            draw_rectangle!(pltvec[i], vertices, color="red", alpha=0.5)
        end
        scatter3d!(pltvec[i], [x_ref[1]], [x_ref[2]], [x_ref[3]], color="green")
    end
   
    anim = @animate for i=1:maxanim
        layout = @layout [a b c d]
        for j in eachindex(hists)
            scatter3d!(pltvec[j], [xvec[j][i]], [yvec[j][i]], [zvec[j][i]],
                xlims=(-2,2),
                ylims=(-2,2),
                zlims=(-8,2),
                camera = (20, 20),
                color = hist_colors[j],
                opacity = feasvec[j][i] ? 0.5 : 0.1
            )
        end
        
        plot(pltvec[1], pltvec[2], pltvec[3], pltvec[4], layout=layout, size=(900, 300))
    end
    display(gif(anim, "hex_anim.gif", fps=10))
end