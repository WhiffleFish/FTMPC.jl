Base.@kwdef struct NonlinearHexModel
    J::SVector{3,Float64} = SA[Jx, Jy, Jz]
    m::Float64 = m
end

function f!(du, u, pa, t)
    Jx, Jy, Jz, m, g, ufunc = pa
    T, u_ϕ, u_θ, u_ψ = ufunc
    x, y, z, ϕ, θ, ψ, uu, v, w, p, q, r = u
    du[1] = cos(θ)*cos(ψ)*uu + (cos(ψ)*sin(ϕ)*sin(θ) - cos(ϕ)*sin(ψ))*v + (sin(ϕ)*sin(ψ) + cos(ϕ)*cos(ψ)*sin(θ))*w #x
    du[2] = cos(θ)*sin(ψ)*uu + (cos(ϕ)*cos(ψ) + sin(ϕ)*sin(θ)*sin(ψ))*v + (cos(ϕ)*sin(θ)*sin(ψ) - sin(ϕ)*cos(ψ))*w#y
    du[3] = -sin(θ)*uu + cos(θ)*sin(ϕ)*v + cos(θ)*cos(ϕ)*w#z
    du[4] = p + (sin(ϕ)*tan(θ))*q + (cos(ϕ)*tan(θ))*r#ϕ
    du[5] = cos(ϕ)*q - sin(ϕ)*r#θ
    du[6] = (sin(ϕ)*sec(θ))*q + (cos(ϕ)*sec(θ))*r#ψ
    du[7] = v*r - q*w - g*sin(θ)#u
    du[8] = p*w - r*uu + g*sin(ϕ)*cos(θ) #v
    du[9] = q*uu - p*v + g*cos(ϕ)*cos(θ) - T/m  #w
    du[10] = ((Jy - Jz)/Jx)*q*r + u_ϕ/Jx  #p
    du[11] = ((Jz - Jx)/Jy)*p*r + u_θ/Jy  #q
    du[12] = ((Jx - Jy)/Jz)*p*q + u_ψ/Jz  #r
end
