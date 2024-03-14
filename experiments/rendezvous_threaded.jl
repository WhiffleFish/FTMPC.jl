using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Clarabel
using DataFrames
using JLD2

function setup(nval)
    ground = 10
    side = 10
    γside = 1.0#0.5e-1
    γground = 1.0#0.5e-1

    constraints = MPC.RendezvousBarrier(h=ground, w=side, l=side, γg = γground, γs = γside)

    modes = 0:1
    T = 30
    Δt = 10#0.05

    u_bounds = (-0.1, 0.1)

    num_modes = length(modes)

    nvec = [0.056, nval]
    models = [MPC.RendezvousModel(nvec[i]) for i in 1:num_modes]
    sys = MPC.BatchDynamics(models; T, Δt, u_bounds)

    ns = sys.inner_statedim
    nm = sys.inner_controldim
    x0 = zeros(ns)
    x0[1:3] .= [0, 8, 0]
    x_ref = zeros(ns)
    x_ref[1:3] .= [0, 2, 0]
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
        verbose = false,
        max_iter= 50_000
    )
    model = JuMPModel(f, x0)

    return model, f, x0

end

function run_sim(simtime, failtime, failmode, delaytime, model, f, x0)    
    unit_planner = MPC.FTMPCPlanner(model, f, 1)
    unit_sim = Simulator(unit_planner, x0=x0, T=simtime, failure=MPC.FixedFailure(failtime,failmode;delay=delaytime))
    unit_hist = simulate(unit_sim)
    return unit_hist
    #p = plot(unit_hist, Δt, side, ground)
end

function run_simulations()
    simtime = 40
    failtime = 10
    failtimes = 1:5:simtime
    numfailtimes = length(failtimes)
    failmode = 2
    delaytime = 1
    ndelays = 5
    numnvals = 10
    nvals = 0.056:-0.001:0.056-0.001*numnvals
    histvec = Vector{MPC.ModeChangeSimHist}(undef, numfailtimes*(numnvals+1)*(ndelays+1))
    #histvec = Vector{MPC.ModeChangeSimHist}()
    histcount = 1
    Threads.@threads for nval in nvals
        Threads.@threads for failtime ∈ failtimes
            delaytimes = failtime:failtime + ndelays
            Threads.@threads for delaytime ∈ delaytimes
                model, f, x0 = setup(nval)
                histvec[histcount] = run_sim(simtime, failtime, failmode, delaytime, model, f, x0)
                histcount += 1
            end
        end
    end

    delaytimes = reduce(vcat, collect.([failtime:failtime+2 for failtime ∈ failtimes]))
    return histvec, simtime, nvals, failtimes, ndelays
end

hists, simtime, nvals, failtimes, ndelays = run_simulations()

jldsave(joinpath(@__DIR__,"results/rendezvous_threaded.jld2"), hists=hists, simtime=simtime, nvals=nvals, 
                                            failtimes=failtimes, ndelays=ndelays)

nothing