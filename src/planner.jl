struct FTMPCPlanner{F<:Union{JuMPFormulator, BarrierJuMPFormulator}}
    model::JuMP.Model
    f::F
end

time_step(p::FTMPCPlanner) = time_step(p.f)
n_modes(p::FTMPCPlanner) = n_modes(p.f.sys)
horizon(p::FTMPCPlanner) = horizon(p.f.sys)

function set_initialstate(p::FTMPCPlanner, x)
    set_initialstate(p.model, p.f.sys, x)
end

function action(p::FTMPCPlanner, x::AbstractVector)
    nm, T = n_modes(p), horizon(p)
    set_initialstate(p, x)
    optimize!(p.model)
    nx,nu = size(p.f.sys.B)
    return value.(p.model[:x][nx+1 : nx+HEX_U_DIM])
end

##

struct ConsensusSearchPlanner{F<:Union{JuMPFormulator, BarrierJuMPFormulator}}
    model::JuMP.Model
    f::F
end

function set_consensus_horizon(p::FTMPCPlanner, t::Int) end
