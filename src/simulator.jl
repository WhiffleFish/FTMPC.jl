# https://jump.dev/JuMP.jl/stable/manual/objective/#Modify-an-objective
struct Simulator{P<:FTMPCPlanner}
    sys::LinearHexModel
    planner::P
    x0::Vector{Float64}
    T::Int
    progress::Bool
end

function Simulator(sys::LinearHexModel, planner::FTMPCPlanner; x0=zeros(12), T=50, progress=true)
    return Simulator(sys, planner, x0, T, progress)
end

struct SimHist
    x::Matrix{Float64}
    u::Matrix{Float64}
    t::Vector{Float64}
end

function simulate(sim::Simulator)
    (;sys, T, planner) = sim
    Δt = time_step(planner)
    ss = c2d(sys.ss, Δt)
    x = sim.x0

    x_hist = Vector{Float64}[]
    u_hist = Vector{Float64}[]

    prog = Progress(T; enabled=sim.progress)
    for t ∈ 1:T
        u = action(planner, x)
        push!(x_hist, x)
        push!(u_hist, u)
        δu = u - sys.u
        x = dstep(ss, x, δu)
        next!(prog)
    end
    finish!(prog)

    return SimHist(
        reduce(hcat, x_hist),
        reduce(hcat, u_hist),
        collect(0.0:Δt:(T-1)*Δt)
    )
end

# FIXME: Doesn't account for multiple barrier functions
function max_barrier_violation(sim::Simulator, hist::SimHist)
    f = sim.planner.f
    ub = first(f.barrier.ub)
    return maximum(dot(f.barrier.A[1,1:12],x)-ub for x ∈ eachcol(hist.x))
end

@recipe function plot(sim::SimHist)
    layout := (1, 2)
    @series begin
        subplot := 1
        labels --> permutedims(STATE_LABELS[TRANSLATIONAL_STATES])
        sim.t, trans_states(flip_z(sim.x))'
    end
    @series begin
        subplot := 2
        labels --> permutedims(["u$i" for i ∈ 1:6])
        sim.t, sim.u'
    end
end
