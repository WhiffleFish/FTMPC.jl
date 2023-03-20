# https://jump.dev/JuMP.jl/stable/manual/objective/#Modify-an-objective
struct Simulator{P<:FTMPCPlanner}
    sys::LinearHexModel
    planner::P
    x0::Vector{Float64}
    T::Int
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

    for t ∈ 1:T
        u = action(planner, x)
        push!(x_hist, x)
        push!(u_hist, u)
        δu = u - sys.u
        x = dstep(ss, x, δu)
    end

    return SimHist(
        reduce(hcat, x_hist),
        reduce(hcat, u_hist),
        collect(0.0:Δt:T*Δt)
    )
end
