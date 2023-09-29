begin
    using BarrierFTMPC
    using Pkg
    const MPC = BarrierFTMPC
    using JuMP
    using OSQP
    using LinearAlgebra
    using Plots
    using ControlSystemsBase

    default(grid=false, framestyle=:box, fontfamily="Computer Modern", label="")
end

begin

    modes = 0:1
    T = 10
    Δt = 0.05
    c2dΔt = Δt #0.1

    #u_bounds = (-Inf,Inf)
    u_bounds = (-20.,20.)
    

    #x0 = [250., 50., 10., 8.]
    x0 = [10., 10., -10., -10.]
    xnom = [0., 0., 0., 0.]
    
    u0 = [0., 0.]
    
    nm = length(modes)
    nn = length(x0)
    np = length(u0) 

    # nommats = vtol_linear_nominal()
    # nominal_model = ss(nommats...)
    # sysd = c2d(nominal_model, c2dΔt)

    # faild = deepcopy(sysd)

    # faild.B[1,1] = 0.0045
    # faild.B[2,1] = 0.0
    # faild.B[3,1] = -0.1319
    # faild.B[2,2] = -0.3624
    # faild.B[3,2] = 0.1053

    sysd = MPC.DoubleIntegrator(;ueff=1.0)
    faild = MPC.DoubleIntegrator(;ueff=0.05)

    nommodel = MPC.CTLinearModel(sysd, xnom, u0)
    failmodel = MPC.CTLinearModel(faild, xnom, u0)

    models = [nommodel, failmodel] # discrete models to batch
end
# Get DT linear models for each mode


sys = MPC.BatchDynamics(models, c2dΔt; T, u_bounds)
#constraints = MPC.VTOLCeiling(sys; h=60, γ=1e-1)

h = -0.0; γ = 1e-1
#γ = 1.0

# positive h and basis → ax <= b

constraint1 = LinearConstraint(-basis(MPC.inner_statedim(sys), 2)*1, h, γ)
constraints = [constraint1]
#constraints = []
x_ref = [0., 0., 0., 0.]
#x_ref = [0., 0., 0., 0.]
Q = I(MPC.inner_statedim(sys))
R = I(MPC.inner_controldim(sys))*1e-1
#Q[2,2] = 20.0
f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    Q=Q,
    R=R,
    constraints,
    eps_prim_inf = 1e-3,
    eps_abs = 1e-5,
    eps_rel = 1e-5,
    verbose = false,
    max_iter= 50_000
)
model = JuMPModel(f, x0)
unit_planner = MPC.FTMPCPlanner(model, f, 1)
a, info = MPC.action_info(unit_planner, x0)


failT = 0.1
simT = 350
failT = Int(round(failT*simT))
delayT = 20
unit_sim = Simulator(unit_planner, x0=x0, T=simT, failure=MPC.FixedFailure(failT,2;delay=delayT))
unit_hist = simulate(unit_sim)
display(plot(unit_hist, ["x", "y", "ẋ", "ẏ"], true, failT, delayT, lw=2, title="Unit"))


consensus_planner = MPC.ConsensusSearchPlanner(model, f)
consensus_sim = Simulator(consensus_planner, x0=x0, T=simT, failure=MPC.FixedFailure(failT,2;delay=delayT))
consensus_hist = simulate(consensus_sim)
display(plot(consensus_hist, ["x", "y", "ẋ", "ẏ"], true, failT, delayT, lw=2, title="Consensus"))
plot([unit_hist.x[1,:], consensus_hist.x[1,:]], [unit_hist.x[2,:],consensus_hist.x[2,:]], 
    [unit_hist.t, consensus_hist.t], failT, delayT, lw=3, title="Full vs One-step Trajectory")
#a1, info1 = MPC.action_info(consensus_planner, x0) 
#print("DONE")