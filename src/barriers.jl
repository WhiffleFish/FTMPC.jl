struct LinearConstraint
    A::Vector{Float64}
    b::Float64
    γ::Float64
    LinearConstraint(A,b,γ=1e-2) = new(A,b,γ)
end

struct BarrierConstraint{MAT<:AbstractMatrix, L<:AbstractVector, U<:AbstractVector}
    A::MAT
    lb::L
    ub::U
end

function linear_constraint_barrier(sys,a,b,γ)
    nm, T = n_modes(sys), horizon(sys)

    A = zeros(T*nm, T*nm*12 + (T-1)*nm*6)

    for t ∈ 1:T-1
        t_section = (t-1)*nm*12
        next_t_section = t*nm*12
        for mode ∈ 1:nm
            x_section = t_section + (mode-1)*12
            x_idxs = x_section+1 : x_section+12
            next_x_section = next_t_section + (mode-1)*12
            next_x_idxs = next_x_section+1 : next_x_section+12

            section_i = (t-1)*nm + mode
            @. A[section_i, next_x_idxs] .= -a
            @. A[section_i, x_idxs] = (1-γ)*a
        end
    end
    ub = fill(γ*b, T*nm)
    lb = fill(-Inf, T*nm)
    return BarrierConstraint(A, lb, ub)
end

struct BarrierJuMPFormulator{T1,T2,T3,T4,T5,T6,T7,T8}
    sys::HexBatchDynamics{T1,T2}
    solver::T3
    P::T4
    q::T5
    barrier::T6
    x_ref_full::T7
    kwargs::T8
end

time_step(f::BarrierJuMPFormulator) = time_step(f.sys)

barrier_constraints(sys, c::LinearConstraint) = linear_constraint_barrier(sys,c.A,c.b,c.γ)

function barrier_constraints(sys, cs::AbstractVector{LinearConstraint})
    As = Matrix{Float64}[]
    lbs = Vector{Float64}[]
    ubs = Vector{Float64}[]
    for c ∈ cs
        barrier = linear_constraint_barrier(sys,c.A,c.b,c.γ)
        push!(As, barrier.A)
        push!(lbs, barrier.lb)
        push!(ubs, barrier.ub)
    end
    return BarrierConstraint(
        reduce(vcat, As),
        reduce(vcat, lbs),
        reduce(vcat, ubs)
    )
end

barrier_constraints(sys,::Nothing) = nothing

function BarrierJuMPFormulator(
    sys::HexBatchDynamics,
    solver;
    P=I(12),
    Q=I(6),
    x_ref=zeros(12),
    constraints = nothing,
    kwargs...)
    # @assert size(P)     == (12,12)
    # @assert size(Q)     == (6,6)
    @assert size(x_ref) == (12,)

    nm, T = n_modes(sys), horizon(sys)
    P_full = process_P(P, nm, T)
    Q_full = process_P(Q, nm, T-1)
    x_ref_full = repeat(x_ref, nm*T)

    P_osqp = blkdiag((P_full, Q_full))
    q_osqp = vcat(vec(-x_ref_full' * P_full), zeros(size(Q_full,1)))

    barrier = barrier_constraints(sys,constraints)

    return BarrierJuMPFormulator(sys, solver, P_osqp, q_osqp, barrier, x_ref_full, convert_kwargs(kwargs))
end

function JuMPModel(f::BarrierJuMPFormulator, x0)
    (;sys, barrier) = f
    (;A,B,Δ_nom,u_bounds) = sys
    u_lower, u_upper = u_bounds

    nm, T = n_modes(sys), horizon(sys)
    model = Model(
        optimizer_with_attributes(f.solver, f.kwargs...),
    )
    nx, nu = size(B)

    A_eq = [I(nx) -B]
    eq_rhs = A*repeat(x0, nm) - Δ_nom

    C = consensus_constraint(sys)
    consensus_rhs = zeros(size(C,1))

    @variable(model, x[1:nx+nu])
    @constraint(model, DYNAMICS, A_eq*x .== eq_rhs)
    @constraint(model, CONSENSUS, C*x[nx+1:end] .== consensus_rhs)
    if u_bounds ≠ (-Inf, Inf)
        @constraint(model, CONTROL, u_lower .≤ x[nx+1:end] .≤ u_upper)
    end
    if !isnothing(barrier)
        @constraint(model, BARRIER, barrier.A*x .≤ barrier.ub)
    end

    @objective(model, Min, 0.5*dot(x, f.P, x) + dot(f.q,x))

    return model
end

struct HexOSQPResults{T}
    X::Vector{Matrix{Float64}}
    U::Vector{Matrix{Float64}}
    t::T
end

function HexOSQPResults(f::BarrierJuMPFormulator, model::JuMP.Model)
    sys = f.sys
    Δt = sys.Δt
    nm, T = n_modes(sys), horizon(sys)
    X = value.(model[:x])
    x = X[1:HEX_X_DIM*nm*T]
    u = X[HEX_X_DIM*nm*T+1 : end]

    X = unbatch_and_disjoint(x, nm, T, HEX_X_DIM)
    U = unbatch_and_disjoint(u, nm, T-1, HEX_U_DIM)
    t = 0.0:Δt:Δt*(T-1)
    return HexOSQPResults(X,U,t)
end

function max_barrier_violation(f::BarrierJuMPFormulator, model::JuMP.Model)
    x = value.(model[:x])
    return maximum(f.barrier.A*x .- f.barrier.ub)
end

function optimizer_max_barrier_violation(model::JuMP.Model)
    return maximum(
        value.(model[:BARRIER]) .- normalized_rhs.(model[:BARRIER])
    )
end

function set_consensus_horizon(model::JuMP.Model, sys::HexBatchDynamics, T::Int)
    (;A,B,Δ_nom,u_bounds) = sys
    nm, T = n_modes(sys), horizon(sys)
    nx, nu = size(B)
end

function set_objective_weights(model::JuMP.Model, f, ws::AbstractVector)
    nm,T = n_modes(f.sys), horizon(f.sys)
    @assert length(ws) == nm
    x = model[:x]
    
    #FIXME: shouldn't be hardcoding P,Q
    P = [I(12)*w for w ∈ ws] 
    Q = [I(6)*w for w ∈ ws]
    
    P_full = process_P(P, nm, T)
    Q_full = process_P(Q, nm, T-1)
    P_osqp = blkdiag((P_full, Q_full))
    q_osqp = vcat(vec(-f.x_ref_full' * P_full), zeros(size(Q_full,1)))

    @objective(model, Min, 0.5*dot(x, P_osqp, x) + dot(q_osqp,x))
    return model
end

function set_initialstate(model::JuMP.Model, sys::HexBatchDynamics, x0)
    (;A,Δ_nom) = sys
    nm = n_modes(sys)
    set_normalized_rhs.(model[:DYNAMICS], A*repeat(x0, nm) - Δ_nom)
    return model
end

function consensus_constraint(sys::HexBatchDynamics, T_consensus=horizon(sys)-1)
    nm, T = n_modes(sys), horizon(sys)
    @assert 1 ≤ T_consensus ≤ T-1
    nu = size(sys.B, 2)

    m = zeros(nm*T_consensus*6, nu)
    for t ∈ 1:T_consensus
        t_section = (t-1)*nm*6
        for mode ∈ 2:nm
            mode_section = t_section + (mode-1)*6
            for u_i ∈ 1:6
                u_idx = mode_section + u_i
                m[u_idx, u_idx] = -1
                m[u_idx, u_idx-6] = 1
            end
        end
    end
    return m
end

function process_P(t::Tuple{<:AbstractMatrix, Int}, nm, T)
    m, idx = t
    ms = fill(zero(m), nm)
    ms[idx] = m
    return process_P(ms, nm, T)
end

function process_P(P::AbstractMatrix, nm, T)
    return sparse(blkdiag(P, nm*T))
end

function process_P(ms::Vector{<:AbstractMatrix}, nm, T)
    return if length(ms) == nm
        sparse(blkdiag(blkdiag(ms),T))
    elseif length(ms) == nm*T
        sparse(blkdiag(ms))
    else
        error("invalid size")
    end
end
