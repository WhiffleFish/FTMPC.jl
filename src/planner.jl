abstract type SimplePlanner end
abstract type ConsensusPlanner end

struct FTMPCPlanner{F<:BarrierJuMPFormulator} <:SimplePlanner
    model::JuMP.Model
    f::F
end

mutable struct NonRobustPlanner{F<:BarrierJuMPFormulator} <:SimplePlanner
    model::JuMP.Model
    f::F
end

mutable struct UnitaryConsensusPlanner{F<:BarrierJuMPFormulator} <:ConsensusPlanner
    model::JuMP.Model
    f::F
end


time_step(p::FTMPCPlanner) = time_step(p.f)
n_modes(p::FTMPCPlanner) = n_modes(p.f.sys)
horizon(p::FTMPCPlanner) = horizon(p.f.sys)


function set_initialstate(p::SimplePlanner, x)
    set_initialstate(p.model, p.f.sys, x)
end

#= function set_initialstate(p::NonRobustPlanner, x)
    set_initialstate(p.model, p.f.sys, x)
end =#

function set_objective_weights(p::FTMPCPlanner, ws::AbstractVector)
    set_objective_weights(p.model, p.f.sys, ws)
end

function action(p::FTMPCPlanner, x::AbstractVector)
    set_initialstate(p, x)
    optimize!(p.model)
    nx,nu = size(p.f.sys.B)
    u = value.(p.model[:x][nx+1 : nx+HEX_U_DIM])
    @assert !any(isnan, u)
    return u
end

##

mutable struct ConsensusSearchPlanner{F<:BarrierJuMPFormulator} <:ConsensusPlanner
    model::JuMP.Model
    f::F
end

time_step(p::ConsensusPlanner) = time_step(p.f)
n_modes(p::ConsensusPlanner) = n_modes(p.f.sys)
horizon(p::ConsensusPlanner) = horizon(p.f.sys)

function set_initialstate(p::ConsensusPlanner, x)
    set_initialstate(p.model, p.f.sys, x)
end

set_consensus_horizon(p::ConsensusPlanner, t) = set_consensus_horizon(p.model, p.f, t)

function set_consensus_horizon(model::JuMP.Model, f, t::Int)
    (;sys) = f
    nm, T = n_modes(sys), horizon(sys)
    nx, nu = size(sys.B)
    C = model[:CONSENSUS]
    
    @assert t < T
    @assert length(C) == nu == 6*nm*(T-1)

    u = model[:x][nx+1:end]
    for t ∈ 1:t
        t_section = (t-1)*nm*6
        for mode ∈ 2:nm
            mode_section = t_section + (mode-1)*6
            for u_i ∈ 1:6
                u_idx = mode_section + u_i
                JuMP.set_normalized_coefficient(C[u_idx], u[u_idx], -1)
                JuMP.set_normalized_coefficient(C[u_idx], u[u_idx-6], 1)
            end
        end
    end
    for t ∈ t+1 : T-1
        t_section = (t-1)*nm*6
        for mode ∈ 2:nm
            mode_section = t_section + (mode-1)*6
            for u_i ∈ 1:6
                u_idx = mode_section + u_i
                JuMP.set_normalized_coefficient(C[u_idx], u[u_idx], 0)
                JuMP.set_normalized_coefficient(C[u_idx], u[u_idx-6], 0)
            end
        end
    end
    return model
end

function optimizer_action(model, f)
    nx = size(f.sys.B, 1)
    return value.(model[:x][nx+1 : nx+HEX_U_DIM])
end

function action(p::UnitaryConsensusPlanner, x::AbstractVector)
    (;model, f) = p
    
    set_initialstate(p, x)

    t = 1
    set_consensus_horizon(model, f, t)
    optimize!(model)
    u = optimizer_action(model, f)
    res = (model,u)
    
    if isnothing(res)
        @warn("No feasible action")
        return nothing
    else
        (_,u) = res
        #@assert !any(isnan, u)
        return u,t
    end
end

function action(p::ConsensusSearchPlanner, x::AbstractVector)
    (;model, f) = p
    T = horizon(p)
    
    set_initialstate(p, x)
    s = BinaryConsensusSearch(model, f)

    res, t = binary_search_max(s, valid_consensus, T-1)
    
    if isnothing(res)
        @warn("No feasible action")
        return nothing
    else
        (_,u) = res
        @assert !any(isnan, u)
        return u,t
    end
end

function action(p::NonRobustPlanner, x::AbstractVector)
    (;model, f) = p
    
    set_initialstate(p, x)
    optimize!(model)
    nx,_ = size(f.sys.B)
    
    u = value.(model[:x][nx+1 : nx+HEX_U_DIM])
    #@assert !any(isnan, u)
    return u
end

function modify_objective_weights(p::ConsensusPlanner, ws)
    model = modified_objective_model(p.f, ws)
    p.model = model
    return p
end

function modify_objective(p::SimplePlanner, mode)
    model = modified_objective_model(p.f, mode)
    p.model = model
    return p
end

set_objective_weights(p::ConsensusSearchPlanner, w) = set_objective_weights(p.model, p.f, w)
