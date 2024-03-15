using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using BarrierFTMPC
const MPC = BarrierFTMPC
using JuMP
using OSQP
using LinearAlgebra
using Clarabel
using DataFrames
using JLD2

function setup(x0)
    ground = 6
    side = 1

    γside = 0.5e-1
    γground = 0.5e-1

    constraints = MPC.ElevatorShaft(h=ground, w=side, l=side, γg = γground, γs = γside)

    modes = 0:2
    T = 10
    Δt = 0.05#0.05

    u_bounds = (0.1, 20.0)

    nm = length(modes)

    models = [MPC.HexCTLinearModel(mode) for mode in modes]
    sys = MPC.BatchDynamics(models; T, Δt, u_bounds)

    #x0 = zeros(12)
    x_ref = zeros(12)
    x_ref[1:3] .= [-0.7, 0.7, 5]

    Qcustom = I(12) * 1.0e-1
    Qcustom[1,1] = 50
    Qcustom[2,2] = 50
    Qcustom[3,3] = 50

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

    return model, f, x0

end

function run_sim(simtime, failtime, failmode, delaytime, model, f, x0; planner=:unit)    
    if planner == :unit
        planner = MPC.FTMPCPlanner(model, f, 1)
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
    Δfailtime = 10
    failtimes = 1:Δfailtime:simtime
    numfailtimes = length(failtimes)
    failmode = 2
    ndelays = 2
    delaytimes = 0:ndelays

    x0 = zeros(12)
    pos_sparsity = 3
    posvec = range(-0.7, stop=0.7, length=pos_sparsity)
    x0vec = []

    for i in posvec
        for j in posvec
            new_x0 = copy(x0)
            new_x0[1] = i
            new_x0[2] = j
            push!(x0vec, new_x0)
        end
    end
    pos2d = [x[1:2] for x in x0vec]

    histvec = Vector{MPC.ModeChangeSimHist}(undef, numfailtimes*length(x0vec)*(ndelays+1))
    #histvec = Vector{MPC.ModeChangeSimHist}()
    histcount = 1
    Threads.@threads for x0 in x0vec
        Threads.@threads for failtime ∈ failtimes
            Threads.@threads for delaytime ∈ delaytimes
                model, f, x0 = setup(x0)
                histvec[histcount] = run_sim(simtime, failtime, failmode, delaytime, model, f, x0; planner=planner_type)
                if size(histvec[histcount].x,2) < simtime
                    println("n: $(x0[1:2]), failtime: $failtime, delaytime: $delaytime, histcount: $histcount,  - failed")
                else
                    println("n: $(x0[1:2]), failtime: $failtime, delaytime: $delaytime, histcount: $histcount,  - succeeded")
                end
                histcount += 1
            end
        end
    end

    return histvec, simtime, pos2d, failtimes, ndelays
end

hists, simtime, pos2d, failtimes, ndelays = run_simulations(planner_type=:unit)
hists_con, _, _, _, _ = run_simulations(planner_type=:consensus)


jldsave(joinpath(@__DIR__,"results/hex_threaded_unit.jld2"), hists=hists, simtime=simtime, pos2d=pos2d, 
                                            failtimes=failtimes, ndelays=ndelays)

jldsave(joinpath(@__DIR__,"results/hex_threaded_consensus.jld2"), hists=hists_con, simtime=simtime, pos2d=pos2d, 
                                            failtimes=failtimes, ndelays=ndelays)

nothing