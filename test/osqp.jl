failures = [0,1]
T = 100
Δt = 0.1
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt)
x0 = zeros(12)
x_ref = zeros(12)
x_ref[3] = -5

f = OSQPFormulator(sys;x_ref,Q=I(6)*1e-2, polish=true)
model = OSQPModel(f, x0)
res = OSQP.solve!(model)

x = res.x[1 : T*nm*12]
u = res.x[T*nm*12 + 1 : end]
_x = sys.A*repeat(x0, nm) + sys.B*u - sys.Δ_nom

@test _x ≈ x atol=1e-3
