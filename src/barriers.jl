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
    _nx, _nu = inner_statedim(sys), inner_controldim(sys)

    A = zeros(T*nm, T*nm*_nx + (T-1)*nm*_nu)

    for t ∈ 1:T-1
        t_section = (t-1)*nm*_nx
        next_t_section = t*nm*_nx
        for mode ∈ 1:nm
            x_section = t_section + (mode-1)*_nx
            x_idxs = x_section+1 : x_section+_nx
            next_x_section = next_t_section + (mode-1)*_nx
            next_x_idxs = next_x_section+1 : next_x_section+_nx

            section_i = (t-1)*nm + mode
            @. A[section_i, next_x_idxs] .= -a
            @. A[section_i, x_idxs] = (1-γ)*a
        end
    end
    
    ub = fill(γ*b, T*nm)
    lb = fill(-Inf, T*nm)
    return BarrierConstraint(A, lb, ub)
end

struct BarrierJuMPFormulator{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    sys::BatchDynamics{T1,T2}
    solver::T3
    P::T4
    q::T5
    P_vec::T6
    Q_vec::T7
    barrier::T8
    x_ref_full::T9
    kwargs::T10
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
    sys::BatchDynamics,
    solver;
    Q=I(inner_statedim(sys)),
    R=I(inner_controldim(sys)),
    x_ref=zeros(inner_statedim(sys)),
    constraints = nothing,
    kwargs...)
    @assert size(x_ref) == (inner_statedim(sys),)
    
    nm, T = n_modes(sys), horizon(sys)

    Q_vec = cost_vec(Q, nm)
    R_vec = cost_vec(R, nm)
    
    Q_full = process_P(Q, nm, T)
    R_full = process_P(R, nm, T-1)
    x_ref_full = repeat(x_ref, nm*T)

    P_osqp = blkdiag((Q_full, R_full))
    q_osqp = vcat(vec(-x_ref_full' * Q_full), zeros(size(R_full,1)))

    barrier = barrier_constraints(sys,constraints)

    return BarrierJuMPFormulator(sys, solver, P_osqp, q_osqp, Q_vec, R_vec, barrier, x_ref_full, convert_kwargs(kwargs))
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

    C = sparse(consensus_constraint(sys, T-1))
    consensus_rhs = zeros(size(C,1))

    @variable(model, x[1:nx+nu])
    @constraint(model, DYNAMICS, sparse(A_eq)*x .== sparse(eq_rhs))

    #@constraint(model, FLOOR, x[3:12:nx] .≤ 6.0)
    #@constraint(model, SIDES_X, -5.0 .≤ x[1:12:nx] .≤ 5.0)
    #@constraint(model, SIDES_Y, -5.0 .≤ x[2:12:nx] .≤ 5.0)
    #@warn "typeof C: " typeof(C)
    #@warn "typeof Abarrier: " typeof(barrier.A)
    @constraint(model, CONSENSUS, C*x[nx+1:end] .== consensus_rhs)
    if u_bounds ≠ (-Inf, Inf)
        @constraint(model, CONTROL, u_lower .≤ x[nx+1:end] .≤ u_upper)
    end
    if !isnothing(barrier)
        @constraint(model, BARRIER, sparse(barrier.A)*x .≤ barrier.ub)
    end

    @objective(model, Min, 0.5*dot(x, f.P, x) + dot(f.q,x))
    #@objective(model, Min, 0.0*dot(x, f.P, x) + 0.0*dot(f.q,x))

    # Note: Model must be optimized for FULL consensus horizon before setting lower
    # consensus horizons. Otherwise OSQP gets angry about changing sparsity pattens.
    optimize!(model)

    return model
end

struct OSQPResults{T}
    X::Vector{Matrix{Float64}}
    U::Vector{Matrix{Float64}}
    t::T
end

function OSQPResults(f::BarrierJuMPFormulator, model::JuMP.Model)
    sys = f.sys
    Δt = sys.Δt
    nm, T = n_modes(sys), horizon(sys)
    X = value.(model[:x])
    x = X[1:sys.inner_statedim*nm*T]
    u = X[sys.inner_statedim*nm*T+1 : end]

    X = unbatch_and_disjoint(x, nm, T, sys.inner_statedim)
    U = unbatch_and_disjoint(u, nm, T-1, sys.inner_controldim)
    t = 0.0:Δt:Δt*(T-1)
    return OSQPResults(X,U,t)
end

to_vec(v::AbstractVector) = v
to_vec(v) = [v]

Base.getindex(res::OSQPResults, i) = OSQPResults(to_vec(res.X[i]),to_vec(res.U[i]),res.t)

@recipe function plot(res::OSQPResults)
    N = length(res.X)
    @assert N == length(res.U)
    layout := (N, 2)
    for i ∈ 1:N
        @series begin
            subplot := 2*(i-1) + 1#(i,1)
            # labels --> permutedims(STATE_LABELS[TRANSLATIONAL_STATES])
            res.t, pos_states(flip_z(res.X[i]))'
        end
        @series begin
            subplot := 2*(i-1) + 2
            # labels --> permutedims(["u$i" for i ∈ 1:6])
            res.t[1:end-1], res.U[i]'
        end
    end
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

function modified_objective_model(f, ws::AbstractVector)
    (;sys, barrier, P_vec, Q_vec, x_ref_full) = f
    (;A,B,Δ_nom,u_bounds) = sys
    _nx = inner_statedim(sys)
    nm, T = n_modes(sys), horizon(sys)
    u_lower, u_upper = u_bounds
    x0 = zeros(_nx) # will be replaced in `action` call anyways

    P_vec_new = P_vec .* ws
    Q_vec_new = Q_vec .* ws
    P_full = process_P(P_vec_new, nm, T)
    Q_full = process_P(Q_vec_new, nm, T-1)

    P = blkdiag((P_full, Q_full))
    q = vcat(vec(-x_ref_full' * P_full), zeros(size(Q_full,1)))
    
    model = Model(
        optimizer_with_attributes(f.solver, f.kwargs...),
    )
    nx, nu = size(B)

    A_eq = [I(nx) -B]
    eq_rhs = A*repeat(x0, nm) - Δ_nom

    C = consensus_constraint(sys, T-1)
    consensus_rhs = zeros(size(C,1))

    @variable(model, x[1:nx+nu])
    @constraint(model, DYNAMICS, sparse(A_eq)*x .== eq_rhs)
    @constraint(model, CONSENSUS, sparse(C)*x[nx+1:end] .== consensus_rhs)
    #@constraint(model, FLOOR, x[3:12:nx] .≤ 6.0)
    #@constraint(model, SIDES_X, -5.0 .≤ x[1:12:nx] .≤ 5.0)
    #@constraint(model, SIDES_Y, -5.0 .≤ x[2:12:nx] .≤ 5.0)

    if u_bounds ≠ (-Inf, Inf)
        @constraint(model, CONTROL, u_lower .≤ x[nx+1:end] .≤ u_upper)
    end
    if !isnothing(barrier)
        @constraint(model, BARRIER, sparse(barrier.A)*x .≤ barrier.ub)
    end

    @objective(model, Min, 0.5*dot(x, P, x) + dot(q,x))
    #@objective(model, Min, 0.0*dot(x, P, x) + 0.0*dot(q,x))
    
    # Note: Model must be optimized for FULL consensus horizon before setting lower
    # consensus horizons. Otherwise OSQP gets angry about changing sparsity pattens.
    optimize!(model)

    return model
end

function set_objective_weights(model::JuMP.Model, f::BarrierJuMPFormulator, ws::AbstractVector)
    (;sys) = f
    nm,T = n_modes(f.sys), horizon(f.sys)
    _nx = inner_statedim(sys)
    _nu = inner_controldim(sys)
    @assert length(ws) == nm
    x = model[:x]
    
    #FIXME: shouldn't be hardcoding P,Q
    P = [I(_nx)*w for w ∈ ws] 
    Q = [I(_nu)*w for w ∈ ws]
    
    P_full = process_P(P, nm, T)
    Q_full = process_P(Q, nm, T-1)
    P_osqp = blkdiag((P_full, Q_full))
    q_osqp = vcat(vec(-f.x_ref_full' * P_full), zeros(size(Q_full,1)))

    @objective(model, Min, 0.5*dot(x, P_osqp, x) + dot(q_osqp,x))
    return model
end

function set_initialstate(model::JuMP.Model, sys::BatchDynamics, x0)
    (;A,Δ_nom) = sys
    nm = n_modes(sys)
    set_normalized_rhs.(model[:DYNAMICS], A*repeat(x0, nm) - Δ_nom)
    return model
end

function consensus_constraint(sys::BatchDynamics, T_consensus=horizon(sys)-1)
    nm, T = n_modes(sys), horizon(sys)
    @assert 1 ≤ T_consensus ≤ T-1
    _nx = inner_statedim(sys)
    _nu = inner_controldim(sys)
    nu = size(sys.B, 2)

    m = zeros(nu, nu)
    for t ∈ 1:T_consensus
        t_section = (t-1)*nm*_nu
        for mode ∈ 2:nm
            mode_section = t_section + (mode-1)*_nu
            for u_i ∈ 1:_nu
                u_idx = mode_section + u_i
                m[u_idx, u_idx] = -1
                m[u_idx, u_idx-_nu] = 1
            end
        end
    end
    return m
end

function cost_vec(t::Tuple{<:AbstractMatrix, Int}, nm)
    m, idx = t
    ms = fill(zero(m), nm)
    ms[idx] = m
    return ms
end

function cost_vec(m::AbstractMatrix, nm)
    return fill(m, nm)
end

function cost_vec(ms::Vector{<:AbstractMatrix}, nm)
    @assert length(ms) == nm
    return ms
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

function consensus_from_hist()
end
