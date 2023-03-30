function binary_search_max(f, isvalid, T)
    best_t = 0
    _min = 0
    _max = T
    best_val = nothing
    while _min < _max
        t = (_min + _max) ÷ 2
        val = f(t)
        if isvalid(val) # consensus horizon too low
            _min = t + 1
            best_t = t
            best_val = val
        else            # consensus horizon too high
            _max = t - 1
        end
    end
    max_val = f(_max)
    return if isvalid(max_val)
        max_val, _max
    else
        best_val, best_t
    end
end

const VALID_STATUSES = [
    OPTIMAL,
    LOCALLY_SOLVED,
    ALMOST_OPTIMAL,
    ALMOST_LOCALLY_SOLVED
]

const WARN_STATUSES = [
    ITERATION_LIMIT,
    TIME_LIMIT,
    NODE_LIMIT,
    SOLUTION_LIMIT,
    MEMORY_LIMIT,
    OBJECTIVE_LIMIT,
    NORM_LIMIT,
    OTHER_LIMIT
]

const INVALID_STATUSES = [
    INFEASIBLE,
    DUAL_INFEASIBLE,
    ALMOST_INFEASIBLE,
    ALMOST_DUAL_INFEASIBLE,
    LOCALLY_INFEASIBLE
]

const ERROR_STATUSES = [
    OPTIMIZE_NOT_CALLED,
    NUMERICAL_ERROR,
    INVALID_MODEL,
    INVALID_OPTION,
    INTERRUPTED,
    OTHER_ERROR
]

struct BinaryConsensusSearch{F<:BarrierJuMPFormulator}
    model::JuMP.Model
    f::F
end

function (b::BinaryConsensusSearch)(t)
    (;model, f) = b
    @show t
    set_consensus_horizon(model, f, t)
    optimize!(model)
    return model, optimizer_action(model, f)
end

function valid_consensus((model,u)::Tuple{JuMP.Model, <:AbstractVector})
    status = termination_status(model)
    if status ∈ VALID_STATUSES
        @assert !any(isnan, u) string(status)
        return true
    elseif status ∈ WARN_STATUSES
        @warn(status)
        return false # FIXME: ehhhhhh?
    elseif status ∈ INVALID_STATUSES
        return false
    else
        error(status)
    end 
end



#=
Enum MathOptInterface.TerminationStatusCode:
OPTIMIZE_NOT_CALLED = 0
OPTIMAL = 1
INFEASIBLE = 2
DUAL_INFEASIBLE = 3
LOCALLY_SOLVED = 4
LOCALLY_INFEASIBLE = 5
INFEASIBLE_OR_UNBOUNDED = 6
ALMOST_OPTIMAL = 7
ALMOST_INFEASIBLE = 8
ALMOST_DUAL_INFEASIBLE = 9
ALMOST_LOCALLY_SOLVED = 10
ITERATION_LIMIT = 11
TIME_LIMIT = 12
NODE_LIMIT = 13
SOLUTION_LIMIT = 14
MEMORY_LIMIT = 15
OBJECTIVE_LIMIT = 16
NORM_LIMIT = 17
OTHER_LIMIT = 18
SLOW_PROGRESS = 19
NUMERICAL_ERROR = 20
INVALID_MODEL = 21
INVALID_OPTION = 22
INTERRUPTED = 23
OTHER_ERROR = 24
=#
