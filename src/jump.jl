struct JuMPFormulator{T1,T2,T3,T4,T5}
    sys::HexBatchDynamics{T1,T2}
    solver::T3
    P::T4
    q::T5
end


function JuMPFormulator(sys::HexBatchDynamics, solver; P=I(12), Q=I(6), x_ref=zeros(12), kwargs...)
    @assert size(P)     == (12,12)
    @assert size(Q)     == (6,6)
    @assert size(x_ref) == (12,)

    nm, T = n_modes(sys), horizon(sys)
    P_full = sparse(blkdiag(P, nm*T))
    Q_full = sparse(blkdiag(Q, nm*(T-1)))
    x_ref_full = repeat(x_ref, nm*T)

    P_osqp = blkdiag((P_full, Q_full))
    q_osqp = vcat(vec(-x_ref_full' * P_full), zeros(size(Q_full,1)))

    return JuMPFormulator(sys, solver, P_osqp, q_osqp)
end

function JuMPModel(f::JuMPFormulator, x0)
    sys = f.sys
    (;A,B,Δ_nom,u_bounds) = sys
    u_lower, u_upper = u_bounds

    nm, T = n_modes(sys), horizon(sys)
    model = Model(f.solver)
    nx, nu = size(B)
    # X = Ā*x_0 + B̄*U - Δ_nom
    # X - B*U = Ā*x_0 - Δ_nom

    A_eq = [
        I(nx) -B;
        zeros(nu,nx) consensus_constraint(sys)
    ]
    eq_rhs = [A*repeat(x0, nm) - Δ_nom; zeros(nu)]

    @variable(model, x[1:nx+nu])
    @constraint(model, A_eq*x .== eq_rhs)
    @constraint(model, u_lower .≤ x[nx+1:end] .≤ u_upper)
    # @objective(model, Min, 0.5*x*f.P*x + f.q*x)
    @objective(model, Min, 0.5*dot(x, f.P, x) + dot(f.q,x))


    return model
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
