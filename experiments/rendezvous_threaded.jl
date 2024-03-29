using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Clarabel
using JLD2
using Pkg

function setup(nval)
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
    Δt = 10

    u_bounds = (-0.1, 0.1)

    num_modes = length(modes)

    nvec = [0.061, nval]
    models = [MPC.RendezvousModel(nvec[i]) for i in 1:num_modes]
    sys = MPC.BatchDynamics(models; T, Δt, u_bounds)

    ns = sys.inner_statedim
    nm = sys.inner_controldim
    x0 = zeros(ns)
    x0[1:3] .= [0, 3.8, 0]
    x_ref = zeros(ns)
    x_ref[1:3] .= [0, 0.5, 0]

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
        verbose = false,
        max_iter= 50_000
    )
    model = JuMPModel(f, x0)


    return model, f, x0, x_ref

end

function run_sim(simtime, failtime, failmode, delaytime, model, f, x0; planner=:unit)    
    if planner == :unit
        planner = MPC.FTMPCPlanner(model, f, 1)
        sim = Simulator(planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
        hist = simulate(sim)
    elseif planner == :max
        planner = MPC.FTMPCPlanner(model, f, f.sys.T-1)
        sim = Simulator(planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
        hist = simulate(sim)
    elseif planner == :nonrobust
        planner = MPC.FTMPCPlanner(model, f, 0)
        sim = Simulator(planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
        hist = simulate(sim)
    elseif planner == :consensus
        planner = MPC.ConsensusSearchPlanner(model, f)
        sim = Simulator(planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
        hist = simulate(sim)
    end
    
    return hist
end

function run_simulations(;planner_type=:unit)
    simtime = 40
    failtimes = 1:5:simtime
    numfailtimes = length(failtimes)
    failmode = 2
    ndelays = 5
    delaytimes = 0:ndelays
    numnvals = 10
    meanmotion = 0.061
    Δn = 0.005
    nvals = meanmotion:-Δn:meanmotion-(numnvals*Δn/2)
    histvec = Vector{MPC.ModeChangeSimHist}(undef, numfailtimes*length(nvals)*(ndelays+1))
    histcount = 1
    Threads.@threads for nval in nvals
        Threads.@threads for failtime ∈ failtimes
            Threads.@threads for delaytime ∈ delaytimes
                model, f, x0, _ = setup(nval)
                histvec[histcount] = run_sim(simtime, failtime, failmode, delaytime, model, f, x0; planner=planner_type)
                histcount += 1
            end
        end
    end

    _, _, _, x_ref = setup(nvals[1])

    return histvec, simtime, nvals, failtimes, ndelays, x_ref
end

hists, simtime, nvals, failtimes, ndelays, x_ref = run_simulations(planner_type=:unit)
hists_max, _, _, _, _ = run_simulations(planner_type=:max)
hists_nonrobust, _, _, _, _ = run_simulations(planner_type=:nonrobust)
hists_con, _, _, _, _ = run_simulations(planner_type=:consensus)


jldsave(joinpath(@__DIR__,"results/rendezvous_threaded_unit.jld2"), hists=hists, simtime=simtime, nvals=nvals, 
                                            failtimes=failtimes, ndelays=ndelays, x_ref=x_ref)
            
jldsave(joinpath(@__DIR__,"results/rendezvous_threaded_max.jld2"), hists=hists_max)

jldsave(joinpath(@__DIR__,"results/rendezvous_threaded_nonrobust.jld2"), hists=hists_nonrobust)

jldsave(joinpath(@__DIR__,"results/rendezvous_threaded_consensus.jld2"), hists=hists_con)

nothing