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
    return batch_dynamics(sys.A, sys.B, T)
end

function batch_dynamics(A::AbstractMatrix,B::AbstractMatrix,T::Int)
    nx, nu = B
    Ā = reduce(vcat, A^t for t ∈ 0:T-1)
    ABts = [zero(B)]
    Abt = B
    for t ∈ 1:T-1
        push!(ABts, Abt)
        Abt = A*Abt
    end

    blks = [_f(ABts, t_i, t_j) for t_i ∈ 1:T, t_j ∈ 1:T-1]
    B̄ = mortar(blks)
    return Ā, B̄
end

function unbatch_states(X::AbstractVector, nx::Int)
    T = length(X) ÷ nx
    return reshape(X, nx, T)
end

function unbatch_and_disjoint(X::AbstractVector, n_modes::Int, T::Int, nx::Int)
    @assert length(X) == nx*n_modes*T
    n_seqs = length(X) ÷ n_modes
    mats = Matrix{Float64}[]
    for mode ∈ 1:n_modes
        mat = zeros(nx, T)
        for t ∈ 1:T
            idx0 = n_modes*(t-1)*nx + (mode-1)*nx + 1
            idxf = idx0 + nx-1
            mat[:,t] .= X[idx0:idxf]
        end
        push!(mats, mat)
    end
    return mats
end

##
make_joint(mats::AbstractVector{<:AbstractMatrix}) = cat(mats..., dims=(1,2))
make_joint(mats::AbstractVector{<:AbstractMatrix}, n) = make_joint(mats)
make_joint(mat::AbstractMatrix, n) = cat(Iterators.repeated(mat, n)..., dims=(1,2))

make_joint_b(mats::AbstractVector{<:AbstractMatrix}) = reduce(vcat, mats)
make_joint_b(mats::AbstractVector{<:AbstractMatrix}, n) = make_joint_b(mats)
make_joint_b(mat::AbstractMatrix, n) = reduce(vcat, Iterators.repeated(mat, n))

make_indep_b(mats::AbstractVector{<:AbstractMatrix}) = cat(mats..., dims=(1,2))
make_indep_b(mats::AbstractVector{<:AbstractMatrix}, n) = make_indep_b(mats)
make_indep_b(mat::AbstractMatrix, n) = cat(Iterators.repeated(mat, n)..., dims=(1,2))

make_joint(vs::AbstractVector{<:AbstractVector}) = mortar(vs)
make_joint(vs::AbstractVector{<:AbstractVector}, n) = make_joint(vs)

default_c(B::AbstractVector{<:AbstractMatrix}, n) = default_c(first(B), n)
default_c(B::AbstractMatrix, n) = I(size(B,1)*n)

default_d(B::AbstractVector{<:AbstractMatrix}, n) = default_d(first(B), n)
default_d(B::AbstractMatrix, n) = zeros(size(B,1)*n, size(B,2))

default_d_joint(B::AbstractVector{<:AbstractMatrix}, n) = default_d(first(B), n)
default_d_joint(B::AbstractMatrix, n) = zeros(size(B,1)*n, size(B,2))

"""
Single control vector u operates over all modes jointly
"""
function joint_dynamics(n::Int, A, B, C=default_c(B,n), D=default_d(B,n))
    return make_joint(A, n), make_joint_b(B, n), C, D
end

"""
each mode has its own control vector u
"""
function independent_dynamics(n::Int, A, B, C=default_c(B,n), D=default_d_joint(B,n))
    return make_joint(A, n), make_indep_b(B, n), C, D
end

##

struct HexBatchDynamics{M1<:AbstractMatrix, M2<:AbstractMatrix}
    A::M1
    B::M2
    Δ_nom::Vector{Float64}
    modes::Vector{Int}
    T::Int
end

n_modes(sys::HexBatchDynamics) = length(sys.modes)
horizon(sys::HexBatchDynamics) = sys.T
#=
Joint THEN batch
=#
function HexBatchDynamics(;T=10, Δt=0.1, failures=0:6)
    n = length(failures)
    As = Matrix{Float64}[]
    Bs = Matrix{Float64}[]
    u_noms = Vector{Float64}[]
    for failure in failures
        model = LinearHexModel(failure)
        push!(u_noms, model.u)
        dsys = c2d(model.ss, Δt)
        push!(As, dsys.A)
        push!(Bs, dsys.B)
    end
    A,B,_,_ = independent_dynamics(n, As, Bs)
    Ā, B̄ = batch_dynamics(A,B,T)
    u_nom_t = reduce(vcat, u_noms)
    U_nom = repeat(u_nom_t, T-1)
    Δ_nom = B̄*U_nom

    return HexBatchDynamics(
        sparse(Ā),
        sparse(B̄),
        convert(Vector{Float64}, Δ_nom), # probably shouldn't need to do this
        convert(Vector{Int}, failures),
        T
    )
end
