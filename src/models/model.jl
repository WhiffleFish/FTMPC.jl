struct CTLinearModel
    ss::StateSpace{Continuous, Float64}
    x::Vector{Float64}
    u::Vector{Float64}
    CTLinearModel(ss,x=zeros(size(ss.A, 1)), u=zeros(size(ss.B, 2))) = new(ss,x,u)
end

statedim(model::CTLinearModel) = size(model.ss.A, 1)
controldim(model::CTLinearModel) = size(model.ss.B, 2)
measdim(model::CTLinearModel) = size(model.ss.C, 1)

struct DTLinearModel
    ss::StateSpace{Discrete{Float64}, Float64}
    x::Vector{Float64}
    u::Vector{Float64}
    Δt::Float64
end

DTLinearModel(sys::CTLinearModel, Δt) = DTLinearModel(c2d(sys.ss,Δt), sys.x, sys.u, Δt)

ControlSystemsBase.c2d(sys::CTLinearModel, Δt) = DTLinearModel(sys, Δt)

statedim(model::DTLinearModel) = size(model.ss.A, 1)
controldim(model::DTLinearModel) = size(model.ss.B, 2)
measdim(model::DTLinearModel) = size(model.ss.C, 1)

struct BatchDynamics{M1<:AbstractMatrix,M2<:AbstractMatrix}
    A::M1
    B::M2
    Δ_nom::Vector{Float64}
    modes::Vector{Int}
    T::Int
    Δt::Float64
    u_bounds::NTuple{2,Vector{Float64}}
    inner_statedim::Int
    inner_controldim::Int
    inner_measdim::Int
    u_noms::Vector{Vector{Float64}}
    models::Vector{StateSpace{Discrete{Float64}, Float64}}
end

convert_u_bound(sys, v::AbstractVector) = v
convert_u_bound(sys, v::Number) = fill(v, controldim(sys))

BatchDynamics(model_vec::AbstractVector{<:CTLinearModel}; Δt, kwargs...) = BatchDynamics(
    map(model -> c2d(model, Δt), model_vec); kwargs...
)

"""
Want to be fed a vector of discrete linear models
"""
function BatchDynamics(model_vec::AbstractVector{<:DTLinearModel}; T=10, u_bounds=(-Inf,Inf))
    @assert length(u_bounds) == 2
    # TODO: check that all models have same dims
    n = length(model_vec)
    sys1 = first(model_vec)
    Δt = sys1.Δt
    
    nx = statedim(sys1)
    nu = controldim(sys1)
    ny = measdim(sys1)

    u_lower = convert_u_bound(sys1, first(u_bounds))
    u_lower = repeat(u_lower, n*(T-1))
    u_upper = convert_u_bound(sys1, last(u_bounds))
    u_upper = repeat(u_upper, n*(T-1))

    As = Matrix{Float64}[]
    Bs = Matrix{Float64}[]
    u_noms = Vector{Float64}[]
    discrete_models = StateSpace{Discrete{Float64}, Float64}[]

    for model in model_vec
        push!(u_noms, model.u)
        dsys = model.ss
        # dsys = c2d(model.ss, Δt)
        push!(As, dsys.A)
        push!(Bs, dsys.B)
        push!(discrete_models, dsys)
    end
    
    A,B,_,_ = independent_dynamics(n, As, Bs)
    Ā, B̄ = batch_dynamics(A,B,T)
    u_nom_t = reduce(vcat, u_noms)
    U_nom = repeat(u_nom_t, T-1)
    Δ_nom = B̄*U_nom

    # seems like unnecessary information
    modes = vec(eachindex(model_vec))
    
    return BatchDynamics(
            sparse(Ā), 
            sparse(B̄), 
            convert(Vector{Float64}, Δ_nom),
            convert(Vector{Int}, modes),
            T, 
            Δt, 
            (u_lower, u_upper), 
            nx, nu, ny,
            u_noms,
            discrete_models
        )
end

time_step(sys::BatchDynamics) = sys.Δt
n_modes(sys::BatchDynamics) = length(sys.modes)
horizon(sys::BatchDynamics) = sys.T

inner_statedim(sys::BatchDynamics) = sys.inner_statedim
inner_controldim(sys::BatchDynamics) = sys.inner_controldim
inner_measdim(sys::BatchDynamics) = sys.inner_measdim
