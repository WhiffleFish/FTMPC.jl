failures = [0,1]
T = 50
Δt = 0.1
# u_bounds = (0.,10.)
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[1:3] .= -5

f = JuMPFormulator(sys, OSQP.Optimizer;x_ref,Q=I(6)*1e-2)
model = JuMPModel(f, x0)
optimize!(model)
res = MPC.HexOSQPResults(f, model)

@test res.X[1][:,end] ≈ x_ref atol = 0.1
@test res.X[2][:,end] ≈ x_ref atol = 0.1
@test res.U[1] ≈ res.U[2] atol = 0.1