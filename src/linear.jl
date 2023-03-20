#SS: [x, y, z, dx, dy, dz]
function hex_linear(failure=0)
    A = zeros(nn,nn)
    Bv = zeros(nn,nmFM)

    A[1,7] = 1
    A[2,8] = 1
    A[3,9] = 1
    A[7,5] = -g
    A[8,4] = g
    A[4,10] = 1
    A[5,11] = 1
    A[6,12] = 1

    Bv[9,1] = -1/m
    Bv[10,2] = 1/Jx
    Bv[11,3] = 1/Jy
    Bv[12,4] = 1/Jz

    MixMat = [1                 1           1                  1             1          1       ;
              l/2                l          l/2               -l/2           -l        -l/2      ;
              l*sqrt(3)/2        0      -l*sqrt(3)/2     -l*sqrt(3)/2     0     l*sqrt(3)/2;
              d/b                -d/b           d/b                 -d/b             d/b          -d/b]

    B = Bv*MixMat
    C = Matrix(I(np))
    D = zeros(np, nm)
    if !iszero(failure)
        B[:,failure] .= 0
    end
    return (A,B,C,D)
end

""" Flip direction of z and dz of state vector s.t. +z points up"""
function flip_z end

function flip_z(x::AbstractVector)
    x_out = copy(x)
    x[3] = -x[3]
    x[9] = -x[9]
    return x
end

function flip_z(X::AbstractMatrix)
    X_out = copy(X)
    X_out[3,:] .= -X_out[3,:]
    X_out[9,:] .= -X_out[9,:]
    return X_out
end

"""Get only translational states (x,y,z,dx,dy,dz) from state"""
function trans_states end

function trans_states(x::AbstractVector)
    @assert length(x) == 12
    return x[TRANSLATIONAL_STATES]
end

function trans_states(X::AbstractMatrix)
    @assert size(X,1) == 12
    return X[TRANSLATIONAL_STATES, :]
end

"""
    X = [x, y, z, ϕ, θ, ψ, dx, dy, dz, dϕ, dθ, dψ]

NOTE: x,y,z in inertial frame; rest in body frame
"""
struct LinearHexModel
    "Linearized perturbation model"
    ss::StateSpace{Continuous, Float64}

    "Linearization state x̄ -> x = x̄ + δx"
    x::Vector{Float64}

    "Linearization control ū -> u = ū + δu"
    u::Vector{Float64}

    LinearHexModel(failure::Int=0) = new(
        ss(hex_linear(failure)...),
        zeros(6),
        hover_control(failure)
    )
end

function dstep(dsys::StateSpace{<:Discrete}, x, δu)
    (;A,B) = dsys
    return A*x + B*δu
end

const HEX_X_DIM = 12
const HEX_U_DIM = 6
