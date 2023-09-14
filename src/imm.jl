struct IMM{O}
    modes::Vector{StateSpace{Discrete{Float64}, Float64}}
    u_noms::Vector{Vector{Float64}}
    T::Matrix{Float64}
    weights::Vector{Float64}
    obs_dist::O
end

statedim(sys::StateSpace) = size(sys.B,1)
controldim(sys::StateSpace) = size(sys.B,2)

default_u_nom(sys::StateSpace) = zeros(controldim(sys))
default_u_noms(sys_vec::AbstractVector{<:StateSpace}) = map(default_u_nom, sys_vec)

struct NothingDist end

function IMM(modes, u_noms=default_u_noms(modes); weights=basis(length(modes),1))
    n = length(modes)
    return IMM(
        modes,
        u_noms, # make default zeros
        Array{Float64}(I(n)),
        weights,
        Returns(NothingDist()) # default unused obs_dist
    )
end

Random.rand(rng::Random.AbstractRNG, ::NothingDist) = nothing

function default_imm(planner::Union{ConsensusSearchPlanner, FTMPCPlanner})
    sys = planner.f.sys # BatchDynamics
    return IMM(sys.models, sys.u_noms)
end

@inline modes(imm::IMM) = imm.modes
@inline weights(imm::IMM) = imm.weights

function hex_mode_transition(failure_prob = 0.05)
    T = Matrix{Float64}(I(7))
    T[1,1] = 1-failure_prob
    T[2:end,1] .= failure_prob / 6
    return T
end

function HexIMM(;Δt=0.1, w=basis(7,1), failure_prob=0.05)
    modes = StateSpace{Discrete{Float64}, Float64}[]
    u_noms = Vector{Float64}[]
    for i ∈ 0:6
        sys = LinearHexModel(i)
        mode = c2d(sys.ss,Δt)
        push!(modes, mode)
        push!(u_noms, sys.u)
    end
    # FIXME: horrible obs dist
    default_cov = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.,1.,1.,100, 100, 100]
    obs_dist = (x) -> MvNormal(x, Diagonal(default_cov))
    return IMM(modes, u_noms, hex_mode_transition(failure_prob), w, obs_dist)
end

function update!(imm::IMM, x, u, y)
    (; weights, modes, u_noms, T, obs_dist) = imm
    tmp = T*weights # predictor
    weights .= tmp 
    for mode in eachindex(modes)
        sys = modes[mode]
        u_nom = u_noms[mode]
        δu = u - u_nom
        xp = dstep(sys, x, δu)
        _pdf = pdf(obs_dist(xp), y) # corrector
        weights[mode] *= _pdf
    end
    return normalize!(weights, 1)
end

update!(imm::IMM{Returns{NothingDist}}, x, u, y) = imm.weights
