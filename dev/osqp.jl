using ControlSystems
using Plots
using LinearAlgebra
using BarrierFTMPC
const MPC = BarrierFTMPC

c2d
model = LinearHexModel(0)
Δt = 0.1
dsys =  c2d(model.ss, Δt)
x0 = zeros(12)
Ā, B̄ = MPC.batch_dynamics(dsys, x0, 0, 20)

Ā*x0
