failures = 0:6
nm = length(failures)
T = 11
Δt = 0.1
sys = MPC.HexBatchDynamics(;failures, T, Δt)
x0 = repeat(zeros(12), length(failures))
δu = repeat(ones(6), length(failures)*(T-1))
u_nom = repeat(reduce(vcat,LinearHexModel(i).u for i ∈ failures), T-1)
u = δu + u_nom

x = sys.A*x0 + sys.B*u - sys.Δ_nom
Xs = MPC.unbatch_and_disjoint(x, nm, T, 12)

t = 0:Δt:1
_u(x,t) = ones(6)
for i ∈ 0:6
    sys = LinearHexModel(i).ss
    y, t, x, uout = lsim(sys, _u, t, x0=zeros(12))
    @test Xs[i+1] ≈ x atol=1e-5
end
