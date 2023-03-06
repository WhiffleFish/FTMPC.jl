using BarrierFTMPC
const MPC = BarrierFTMPC
mod = LinearModel()

x0 = rand(size(mod.A,2))
u = rand(size(mod.B,2))

BarrierFTMPC.step(mod, x0, u)

using ControlSystems

ss
hex = ss(MPC.A, MPC.B, MPC.C, MPC.D)

@edit StateSpace(MPC.A, MPC.B, MPC.C, MPC.D)
