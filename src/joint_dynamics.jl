function _f(Abts, t_i, t_j)
    # which diagonal are we on?
    sub_diag = max((t_i - t_j) + 1, 1)
    return Abts[sub_diag]
end

"""
    `batch_dynamics(sys::StateSpace{<:Discrete}, T::Int)`

Given discrete state space system, give Ā, B̄ such that X = Ā*x_init + B̄*U
"""
function batch_dynamics(sys::StateSpace{<:Discrete}, T::Int)
    (;A,B) = sys
    nx, nu = B
    Ā = reduce(vcat, A^t for t ∈ 0:T-1)
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
make_joint(mats::AbstractVector{<:AbstractMatrix}) = cat(mats..., dims=(1,2))
make_joint(mats::AbstractVector{<:AbstractMatrix}, n) = make_joint(mats)
make_joint(mat::AbstractMatrix, n) = cat(Iterators.repeated(mat, n)..., dims=(1,2))

make_joint_b(mats::AbstractVector{<:AbstractMatrix}) = reduce(vcat, mats)
make_joint_b(mats::AbstractVector{<:AbstractMatrix}, n) = make_joint_b(mats)
make_joint_b(mat::AbstractMatrix, n) = reduce(vcat, Iterators.repeated(mat, n))

make_joint(vs::AbstractVector{<:AbstractVector}) = mortar(vs)
make_joint(vs::AbstractVector{<:AbstractVector}, n) = make_joint(vs)

default_c(B::AbstractVector{<:AbstractMatrix}, n) = default_c(first(B), n)
default_c(B::AbstractMatrix, n) = I(size(B,1)*n)

default_d(B::AbstractVector{<:AbstractMatrix}, n) = default_d(first(B), n)
default_d(B::AbstractMatrix, n) = zeros(size(B,1)*n, size(B,2))

function joint_dynamics(n::Int, A, B, C=default_c(B,n), D=default_d(B,n))
    return make_joint(A, n), make_joint_b(B, n), C, D
end

##

# FIXME: horribly inefficient
function hex_batch_dynamics(T=10, Δt=0.1, failures=0:6)
    n = length(failures)
    As = Matrix{Float64}[]
    Bs = Matrix{Float64}[]
    for failure in failures
        sys = LinearHexModel(failure).ss
        push!(As, sys.A)
        push!(Bs, sys.B)
    end
    dsys = c2d(ss(joint_dynamics(n, As, Bs)...), Δt)
    return batch_dynamics(dsys, T)
end
