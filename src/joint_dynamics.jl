function _f(Abts, t_i, t_j)
    # which diagonal are we on?
    sub_diag = max((t_i - t_j) + 1, 1)
    return Abts[sub_diag]
end

"""
    `batch_dynamics(sys::StateSpace{<:Discrete}, T::Int)`

Given discrete state space system
Given linear system, x_init, failure point, and T, get matrix expression for X, J
"""
function batch_dynamics(sys::StateSpace{<:Discrete}, T::Int)
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

function unbatch_states(X::AbstractVector, nx::Int)
    T = length(X) ÷ nx
    return reshape(X, nx, T)
end

##
make_joint(mats::AbstractVector{<:AbstractMatrix}) = cat(mats, dims=(1,2))
make_joint(mats::AbstractVector{<:AbstractMatrix}, n) = make_joint(mats)
make_joint(mat::AbstractMatrix, n) = cat(Iterators.repeated(mat, n)..., dims=(1,2))

make_joint_b(mats::AbstractVector{<:AbstractMatrix}) = reduce(vcat, mats)
make_joint_b(mats::AbstractVector{<:AbstractMatrix}, n) = make_joint(mats)
make_joint_b(mat::AbstractMatrix, n) = reduce(vcat, Iterators.repeated(mat, n))

make_joint(vs::AbstractVector{<:AbstractVector}) = mortar(vs)
make_joint(vs::AbstractVector{<:AbstractVector}, n) = make_joint(vs)

default_c(B::AbstractVector{<:AbstractMatrix}) = default_c(first(B))
default_c(B::AbstractMatrix) = I(size(B,1))

default_d(B::AbstractVector{<:AbstractMatrix}) = default_d(first(B))
default_d(B::AbstractMatrix) = zero(B)

function joint_dynamics(P::AbstractVector{Vector{Float64}}, A, B, C=default_c(B), D=default_d(B))
    n = length(P)
    return make_joint(P), make_joint(A, n), make_joint_b(B, n)
end
