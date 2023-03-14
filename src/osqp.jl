struct OSQPFormulator{T1,T2,T3,T4,T5}
    sys::HexBatchDynamics{T1,T2}
    P::T3 # BIG
    q::T4
    kwargs::T5
end

function OSQPFormulator(sys::HexBatchDynamics; P=I(12), Q=I(6), x_ref=zeros(12), kwargs...)
    @assert size(P)     == (12,12)
    @assert size(Q)     == (6,6)
    @assert size(x_ref) == (12,)

    nm, T = n_modes(sys), horizon(sys)
    P_full = sparse(blkdiag(P, nm*T))
    Q_full = sparse(blkdiag(Q, nm*(T-1)))
    x_ref_full = repeat(x_ref, nm*T)

    P_osqp = blkdiag((P_full, Q_full))
    q_osqp = vcat(vec(-x_ref_full' * P_full), zeros(size(Q_full,1)))

    return OSQPFormulator(sys, P_osqp, q_osqp, kwargs)
end

function consensus_constraint(sys::HexBatchDynamics, T_consensus=horizon(sys)-1)
    nm, T = n_modes(sys), horizon(sys)
    @assert 1 ≤ T_consensus ≤ T-1
    nu = size(sys.B, 2)
    m = zeros(nu, nu)
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

function OSQPModel(f::OSQPFormulator, x0)
    sys = f.sys
    (;A,B,Δ_nom) = sys
    nm, T = n_modes(sys), horizon(sys)
    model = OSQP.Model()
    nx, nu = size(B)
    # X = Ā*x_0 + B̄*U - Δ_nom
    # X - B*U = Ā*x_0 - Δ_nom

    A_constraint = [
        I(nx) -B;
        zeros(nu,nx) consensus_constraint(sys)
    ]
    l = [A*repeat(x0, nm) - Δ_nom; zeros(nu)]
    u = [A*repeat(x0, nm) - Δ_nom; zeros(nu)]

    @assert size(A_constraint, 1) == length(l) == length(u)
    @assert size(A_constraint, 2) == nx + nu == size(f.P,2) == length(f.q)

    OSQP.setup!(model; P=f.P, q=f.q, A=A_constraint, l=l, u=u, f.kwargs...)
    return model
end


struct HexOSQPResults
    X::Vector{Matrix{Float64}}
    U::Vector{Matrix{Float64}}
end

function HexOSQPResults(f::OSQPFormulator, res::OSQP.Results)
    sys = f.sys
    nm, T = n_modes(sys), horizon(sys)
    x = res.x[1:12*nm*T]
    u = res.x[12*nm*T+1 : end]

    X = unbatch_and_disjoint(x, nm, T, HEX_X_DIM)
    U = unbatch_and_disjoint(u, nm, T-1, HEX_U_DIM)
    return HexOSQPResults(X,U)
end
