# https://jump.dev/JuMP.jl/stable/manual/objective/#Modify-an-objective
struct Simulator{P}
    sys::LinearHexModel
    planner::P
    x0::Vector{Float64}
    T::Int
    progress::Bool
end

function Simulator(sys::LinearHexModel, planner; x0=zeros(12), T=50, progress=true)
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
        isnothing(u) && break
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

##
struct FixedFailure
    t::Int
    mode::Int
    instant_update::Bool
    FixedFailure(t, mode; instant_update=true) = new(t,mode, instant_update)
end

function fail!(imm::IMM, planner, t, mode_idx, failure::FixedFailure)
    if t == failure.t
        if failure.instant_update
            imm.weights .= 0.0
            imm.weights[failure.mode] = 1.0
            modify_objective_weights(planner, imm.weights) 
        end
        return failure.mode
    else
        return mode_idx
    end
end

struct NoFailure end

function fail!(imm, planner, t, mode_idx, failure::NoFailure)
    return mode_idx
end

struct RandFailure end

struct ModeChangeSimulator{P,F}
    imm::IMM
    planner::P
    x0::Vector{Float64}
    T::Int
    failure::F
    progress::Bool
end

const FloatRange = StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}

struct ModeChangeSimHist
    x::Matrix{Float64}
    u::Matrix{Float64}
    t::Vector{Float64}
    mode::Vector{Int}
    w::Matrix{Float64}
    info::Vector{HexOSQPResults{FloatRange}}
    consensus::Vector{Int}
end

function Simulator(imm::IMM, planner; x0=zeros(12), failure=NoFailure(), T=50, progress=true)
    return ModeChangeSimulator(imm, planner, x0, T, failure, progress)
end


function simulate(sim::ModeChangeSimulator)
    (;imm, T, planner, failure) = sim
    mode_idx = weighted_sample(imm.weights)
    ss = imm.modes[mode_idx]
    Δt = time_step(planner)
    x = sim.x0

    x_hist = Vector{Float64}[]
    u_hist = Vector{Float64}[]
    mode_hist = Int[]
    imm_state_hist = Vector{Float64}[]
    info = HexOSQPResults{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}[]
    consensus_hist = Int[]

    prog = Progress(T; enabled=sim.progress)
    for t ∈ 1:T
        # FIXME: Set objective weights in sim loop without changing sparsity pattern
        # modify_objective_weights(planner, imm.weights) 
        push!(x_hist, copy(x))
        u,_info,h_star = if planner isa ConsensusSearchPlanner
            action_info(planner, x)
        elseif planner isa FTMPCPlanner
            action_info(planner, x)..., 1
        end
        if isnothing(u)
            @show x
            break
        end
        push!(consensus_hist, h_star)
        push!(u_hist, copy(u))
        push!(mode_hist, mode_idx)
        push!(imm_state_hist, copy(imm.weights))
        push!(info, _info)

        δu = u - imm.u_noms[mode_idx]
        xp = dstep(ss, x, δu)

        y = rand(imm.obs_dist(xp))
        update!(imm, x, u, y)
        x = xp

        # update mode
        # mode_idx = weighted_sample(imm.T[:,mode_idx])
        # ss = imm.modes[mode_idx]

        mode_idx = fail!(imm, planner, t, mode_idx, failure)
        ss = imm.modes[mode_idx]

        next!(prog)
    end
    finish!(prog)

    return ModeChangeSimHist(
        reduce(hcat, x_hist),
        reduce(hcat, u_hist),
        collect(0.0:Δt:(T-1)*Δt),
        mode_hist,
        reduce(hcat, imm_state_hist),
        info,
        consensus_hist
    )
end


##

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
        sim.t[1:size(sim.x,2)], trans_states(flip_z(sim.x))'
    end
    @series begin
        subplot := 2
        labels --> permutedims(["u$i" for i ∈ 1:6])
        sim.t[1:size(sim.u,2)], sim.u'
    end
end

@recipe function plot(sim::ModeChangeSimHist)
    layout := (1, 2)
    @series begin
        subplot := 1
        labels --> permutedims(STATE_LABELS[TRANSLATIONAL_STATES])
        sim.t[1:size(sim.x,2)], trans_states(flip_z(sim.x))'
    end
    @series begin
        subplot := 2
        labels --> permutedims(["u$i" for i ∈ 1:6])
        sim.t[1:size(sim.u,2)], sim.u'
    end
end
