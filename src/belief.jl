function _f(Abts, t_i, t_j)
    # which diagonal are we on?
    sub_diag = max((t_i - t_j) + 1, 1)
    return Abts[sub_diag]
end

"""
Given linear system, x_init, failure point, and T, get matrix expression for X, J
"""
function batch_dynamics(sys::StateSpace{<:Discrete}, x_init, failure::Int, T::Int)
    (;A,B) = sys
    nx, nu = B
    Ā = reduce(vcat, A^t for t ∈ 0:T-1)
    # Ā = mortar([A^t for t ∈ 0:T-1])
    ABts = [zero(B)]
    Abt = B
    for t ∈ 1:T-1
        push!(ABts, Abt)
        Abt = A*Abt
    end # still need to account for changing B

    blks = [_f(ABts, t_i, t_j) for t_i ∈ 1:T, t_j ∈ 1:T-1]
    B̄ = mortar(blks)
    return Ā, B̄
end
