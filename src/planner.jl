mutable struct FTMPCPlanner{F<:BarrierJuMPFormulator}
    model::JuMP.Model
    f::F
    consensus_horizon::Int
    function FTMPCPlanner(model, f, consensus_horizon=horizon(f.sys)-1)
        set_consensus_horizon(model, f, consensus_horizon)
        return new{typeof(f)}(model, f, consensus_horizon)
    end
end

time_step(p::FTMPCPlanner) = time_step(p.f)
n_modes(p::FTMPCPlanner) = n_modes(p.f.sys)
horizon(p::FTMPCPlanner) = horizon(p.f.sys)

function set_initialstate(p::FTMPCPlanner, x)
    set_initialstate(p.model, p.f.sys, x)
end

function set_objective_weights(p::FTMPCPlanner, ws::AbstractVector)
    set_objective_weights(p.model, p.f, ws)
end

action(p::FTMPCPlanner, x::AbstractVector) = first(action_info(p, x))

function action_info(p::FTMPCPlanner, x::AbstractVector)
    set_initialstate(p, x)
    optimize!(p.model)
    nx,nu = size(p.f.sys.B)
    u = value.(p.model[:x][nx+1 : nx+p.f.sys.inner_controldim])
    return if termination_status(p.model) ∈ INVALID_STATUSES
        #@warn("No feasible action")
        nothing, nothing
    else
        u, OSQPResults(p.f, p.model)
    end
end

function action_info(p::FTMPCPlanner, x::AbstractVector, ws::AbstractVector)
    set_initialstate(p, x)
    optimize!(p.model)
    nx,nu = size(p.f.sys.B)
    
    if p.consensus_horizon == 0
        idx = findmax(ws)[2]
        ctrl_idx = nx+1+((idx-1)*p.f.sys.inner_controldim) : nx+p.f.sys.inner_controldim+((idx-1)*p.f.sys.inner_controldim)
        u = value.(p.model[:x][ctrl_idx])
    else
        u = value.(p.model[:x][nx+1 : nx+p.f.sys.inner_controldim])
    end
    return if termination_status(p.model) ∈ INVALID_STATUSES
        #@warn("No feasible action")
        idx = findmax(ws)[2]
        ctrl_idx = nx+1+((idx-1)*p.f.sys.inner_controldim) : nx+p.f.sys.inner_controldim+((idx-1)*p.f.sys.inner_controldim)
        u = value.(p.model[:x][ctrl_idx])

        u, OSQPResults(p.f, p.model)
    else
        u, OSQPResults(p.f, p.model)
    end
end

##

mutable struct ConsensusSearchPlanner{F<:BarrierJuMPFormulator}
    model::JuMP.Model
    f::F
end

time_step(p::ConsensusSearchPlanner) = time_step(p.f)
n_modes(p::ConsensusSearchPlanner) = n_modes(p.f.sys)
horizon(p::ConsensusSearchPlanner) = horizon(p.f.sys)

function set_initialstate(p::ConsensusSearchPlanner, x)
    set_initialstate(p.model, p.f.sys, x)
end

set_consensus_horizon(p::ConsensusSearchPlanner, t) = set_consensus_horizon(p.model, p.f, t)

function set_consensus_horizon(model::JuMP.Model, f, t::Int)
    (;sys) = f
    nm, T = n_modes(sys), horizon(sys)
    _nx, _nu = inner_statedim(sys), inner_controldim(sys)
    nx, nu = size(sys.B)
    C = model[:CONSENSUS]
    
    @assert t < T
    @assert length(C) == nu == _nu*nm*(T-1)

    u = model[:x][nx+1:end]
    for t ∈ 1:t
        t_section = (t-1)*nm*_nu
        for mode ∈ 2:nm
            mode_section = t_section + (mode-1)*_nu
            for u_i ∈ 1:_nu
                u_idx = mode_section + u_i
                JuMP.set_normalized_coefficient(C[u_idx], u[u_idx], -1)
                JuMP.set_normalized_coefficient(C[u_idx], u[u_idx-_nu], 1)
            end
        end
    end
    for t ∈ t+1 : T-1
        t_section = (t-1)*nm*_nu
        for mode ∈ 2:nm
            mode_section = t_section + (mode-1)*_nu
            for u_i ∈ 1:_nu
                u_idx = mode_section + u_i
                JuMP.set_normalized_coefficient(C[u_idx], u[u_idx], 0)
                JuMP.set_normalized_coefficient(C[u_idx], u[u_idx-_nu], 0)
            end
        end
    end
    return model
end

function optimizer_action(model, f)
    nx = size(f.sys.B, 1)
    return value.(model[:x][nx+1 : nx+f.sys.inner_controldim])
end

function action_info(p::ConsensusSearchPlanner, x::AbstractVector)
    (;model, f) = p
    T = horizon(p)
    
    s = BinaryConsensusSearch(model, f)
    set_initialstate(p, x)
    res, t = binary_search_max(s, valid_consensus, T-1)
    if isnothing(res)
        #@warn("No feasible action")
        return nothing, nothing, nothing
    else
        (m,u,info) = res
        @assert !any(isnan, u)
        return u, info, t
    end
end

action(p::ConsensusSearchPlanner, x::AbstractVector) = first(action_info(p, x))

function modify_objective_weights(p::Union{ConsensusSearchPlanner, FTMPCPlanner}, ws)
    model = modified_objective_model(p.f, ws)
    p.model = model
    return p
end

set_objective_weights(p::ConsensusSearchPlanner, w) = set_objective_weights(p.model, p.f, w)
