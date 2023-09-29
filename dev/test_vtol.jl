begin
    using BarrierFTMPC
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
    c2dΔt = 0.1

    #u_bounds = (-Inf,Inf)
    u_bounds = (-200.,200.)
    

    #x0 = [250., 50., 10., 8.]
    x0 = [100., 0., 0., 0.]
    xnom = [135., 0., 0., 0.]
    
    u0 = [100., 100.]
    
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


    sysd = MPC.fault_free_vtol()
    # faild = MPC.system_fault_vtol()
    faild = MPC.actuator_fault_vtol()

    nommodel = MPC.DTLinearModel(sysd, xnom, u0, c2dΔt)
    failmodel = MPC.DTLinearModel(faild, xnom, u0, c2dΔt)

    models = [nommodel, failmodel] # discrete models to batch
end
# Get DT linear models for each mode


sys = MPC.BatchDynamics(models; T, u_bounds)
#constraints = MPC.VTOLCeiling(sys; h=60, γ=1e-1)

h = 0; γ = 1e-1
γ = 1.0

# positive h and basis → ax <= b

constraints = LinearConstraint(-basis(MPC.inner_statedim(sys), 4)*1, h, γ)
#constraints = []
x_ref = [150., 0., 0., 0.]
#x_ref = [0., 0., 0., 0.]
Q = I(MPC.inner_statedim(sys))*1.0
Q[1,1] = 20.0
f = BarrierJuMPFormulator(
    sys,
    OSQP.Optimizer;
    x_ref,
    Q=I(MPC.inner_statedim(sys)),
    R=I(MPC.inner_controldim(sys))*1e-3,
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

failT = 0.2
simT = 100
failT = Int(round(failT*simT))
delayT = 15
unit_sim = Simulator(unit_planner, x0=x0, T=simT, failure=MPC.FixedFailure(failT,2;delay=delayT))
unit_hist = simulate(unit_sim)
display(plot(unit_hist, ["Vh", "Vv", "q", "θ"], true, failT, delayT, lw=2, title="Unit"))


consensus_planner = MPC.ConsensusSearchPlanner(model, f)
consensus_sim = Simulator(consensus_planner, x0=x0, T=simT, failure=MPC.FixedFailure(failT,2;delay=delayT))
consensus_hist = simulate(consensus_sim)
plot(consensus_hist, ["Vh", "Vv", "q", "θ"], true, failT, delayT, lw=2, title="Consensus")
