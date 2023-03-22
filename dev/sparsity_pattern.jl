using BarrierFTMPC
const MPC = BarrierFTMPC
using SparseArrays
using JuMP
using OSQP
using LinearAlgebra
using Plots

failures = 0:6
T = 20
Δt = 0.1
u_bounds = (-Inf,Inf)
nm = length(failures)
sys = MPC.HexBatchDynamics(;failures, T, Δt, u_bounds)
x0 = zeros(12)
x_ref = basis(12,3)

# THIS FAILS
f = BarrierJuMPFormulator(sys,OSQP.Optimizer;x_ref)
model = JuMPModel(f, x0)
MPC.set_consensus_horizon(model, f, 10)
optimize!(model)
MPC.set_consensus_horizon(model, f, 15)

# THIS ALSO FAILS
model = JuMPModel(f, x0)
MPC.set_consensus_horizon(model, f, 15)
optimize!(model)
MPC.set_consensus_horizon(model, f, 10)

# THIS WORKS?!?
model = JuMPModel(f, x0)
optimize!(model)
MPC.set_consensus_horizon(model, f, 15)
optimize!(model)
MPC.set_consensus_horizon(model, f, 10)
optimize!(model)

#= 
NOTE: changing sparsity pattern issue is not raised if model is optimized
for full consensus horizon.
=#
