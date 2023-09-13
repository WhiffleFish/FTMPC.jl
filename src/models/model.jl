abstract type Model end

function statedim(::Model) end

function controldim(::Model) end

function measdim(::Model) end

struct CTLinearModel{AT<:AbstractMatrix, BT<:AbstractMatrix, CT<:AbstractMatrix, DT<:AbstractMatrix} <: Model
    ss::StateSpace{Continuous, Float64}
    x::Vector{Float64}
    u::Vector{Float64}
    CTLinearModel(A::AT,B::BT,C::CT=default_c(B), D::DT=default_d(B)) where {AT,BT,CT,DT} = new{AT,BT,CT,DT}(A,B,C,D)
end

statedim(model::CTLinearModel) = size(model.ss.A, 1)
controldim(model::CTLinearModel) = size(model.ss.B, 2)
measdim(model::CTLinearModel) = size(model.ss.C, 1)

struct BatchDynamics{M1,M2} <: Model
    A::M1
    B::M2
    Δ_nom::Vector{Float64}
    modes::Vector{Int}
    T::Int
    Δt::Float64
    u_bound::Vector{Float64}
    inner_statedim::Int
    inner_controldim::Int
    inner_measdim::Int
end

convert_u_bound(sys, v::AbstractVector) = v
convert_u_bound(sys, v::Number) = fill(v, controldim(sys))

"""
Want to be fed a vector of discrete linear models
"""
function BatchDynamics(model_vec::AbstractVector{<:CTLinearModel}; T=10, Δt=0.1, u_bound=(-Inf,Inf))
    @assert length(u_bound) == 2
    # TODO: check that all models have same dims
    n = length(model_vec)
    sys1 = first(model_vec)
    
    nx = statedim(sys1)
    nu = controldim(sys1)
    ny = measdim(sys1)

    u_lower = convert_u_bound(sys1, first(u_bound))
    u_upper = convert_u_bound(sys1, last(u_bound))
    
    As = Matrix{Float64}[]
    Bs = Matrix{Float64}[]
    u_noms = Vector{Float64}[]

    for model in model_vec
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

    return BatchDynamics(A, B, Δ_nom, modes, T, Δt, (u_lower, u_upper), nx , nu, ny)
end
