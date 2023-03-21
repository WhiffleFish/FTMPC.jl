struct JuMPFormulator{T1,T2,T3,T4,T5,T6}
    sys::HexBatchDynamics{T1,T2}
    solver::T3
    P::T4
    q::T5
    kwargs::T6
end

time_step(f::JuMPFormulator) = time_step(f.sys)

convert_kwargs(kwargs::Base.Pairs) = Tuple(string(a)=>b for (a,b) ∈ kwargs)

function JuMPFormulator(sys::HexBatchDynamics, solver; P=I(12), Q=I(6), x_ref=zeros(12), kwargs...)
    @assert size(P)     == (12,12)
    @assert size(Q)     == (6,6)
    @assert size(x_ref) == (12,)

    nm, T = n_modes(sys), horizon(sys)
    P_full = process_P(P, nm, T)
    Q_full = process_P(Q, nm, T-1)
    x_ref_full = repeat(x_ref, nm*T)

    P_osqp = blkdiag((P_full, Q_full))
    q_osqp = vcat(vec(-x_ref_full' * P_full), zeros(size(Q_full,1)))

    return JuMPFormulator(sys, solver, P_osqp, q_osqp, convert_kwargs(kwargs))
end

function JuMPModel(f::JuMPFormulator, x0)
    sys = f.sys
    (;A,B,Δ_nom,u_bounds) = sys
    u_lower, u_upper = u_bounds

    nm, T = n_modes(sys), horizon(sys)
    model = Model(
        optimizer_with_attributes(f.solver, f.kwargs...),
    )
    nx, nu = size(B)
    # X = Ā*x_0 + B̄*U - Δ_nom
    # X - B*U = Ā*x_0 - Δ_nom

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
    @objective(model, Min, 0.5*dot(x, f.P, x) + dot(f.q,x))

    return model
end

function set_initialstate(model::JuMP.Model, sys::HexBatchDynamics, x0)
    (;A,Δ_nom) = sys
    nm = n_modes(sys)
    set_normalized_rhs.(model[:DYNAMICS], A*repeat(x0, nm) - Δ_nom)
    return model
end


function set_consensus_horizon(model::JuMP.Model, sys::HexBatchDynamics, T::Int)
    (;A,B,Δ_nom,u_bounds) = sys
    nm, T = n_modes(sys), horizon(sys)
    nx, nu = size(B)
end

function HexOSQPResults(f::JuMPFormulator, model::JuMP.Model)
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
