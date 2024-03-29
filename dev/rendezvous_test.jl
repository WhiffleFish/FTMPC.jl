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
    γside = 0.5e-1
    γground = 0.5e-1

    constraints = MPC.RendezvousBarrier(h=ground, w=side, l_lower=ylower, l_upper=yupper,
                                          γg = γground, γs = γside)

    modes = 0:1
    T = 30
    Δt = 7

    u_bounds = (-0.1, 0.1)

    num_modes = length(modes)

    nvec = [0.061, 0.09]
    
    #nvec = [0.061, 0.101]
    models = [MPC.RendezvousModel(nvec[i]) for i in 1:num_modes]
    sys = MPC.BatchDynamics(models; T, Δt, u_bounds)

    ns = sys.inner_statedim
    nm = sys.inner_controldim
    x0 = zeros(ns)
    x0[1:3] .= [0.01, 3.8, 0]
    x_ref = zeros(ns)
    x_ref[1:3] .= [1.0, 1.0, 0]

    Qcustom = I(ns) * 1.0e-1
    Qcustom[1,1] = 50
    Qcustom[2,2] = 50

    Qcustom[3,3] = 50

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
        #iterative_refinement_abstol = 1e-9,
        #iterative_refinement_reltol = 1e-8,
        verbose = false,
        max_iter= 50_000
    )
    model = JuMPModel(f, x0)

    simtime = 40
    failtime = 5
    failmode = 2
    delaytime = 1
end

begin "Unit"
    unit_planner = MPC.FTMPCPlanner(model, f, 1)
    #a, info = MPC.action_info(unit_planner, rand(12)*1e-2)

    unit_sim = Simulator(unit_planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
    unit_hist = simulate(unit_sim)
    #plot(unit_hist, lw=2)
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
    #consensus_hist = simcon(consensus_sim)
    #plot(consensus_hist, lw=2, true)
    #p = plot(consensus_hist, Δt, side, ground)
    barriers = [ylower, yupper]
    params = [Δt, barriers, x_ref, Int(failtime), Int(delaytime)]
    hists = [consensus_hist, unit_hist, nonrobust_hist, full_hist]
    totalplt = plot(hists, params) |> display
    #plot(consensus_hist.consensus) |> display
end
