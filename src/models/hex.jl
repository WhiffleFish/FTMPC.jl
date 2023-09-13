include("hex_constants.jl")

struct HexBatchDynamics{M1<:AbstractMatrix, M2<:AbstractMatrix}
    A::M1
    B::M2
    Δ_nom::Vector{Float64}
    modes::Vector{Int}
    T::Int
    Δt::Float64
    u_bounds::Tuple{Float64, Float64}
end

time_step(sys::HexBatchDynamics) = sys.Δt
n_modes(sys::HexBatchDynamics) = length(sys.modes)
horizon(sys::HexBatchDynamics) = sys.T
#=
Joint THEN batch
=#
function HexBatchDynamics(;T=10, Δt=0.1, failures=0:6, u_bounds=(-Inf,Inf))
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
        T,
        Δt,
        u_bounds
    )
end
