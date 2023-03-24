using Distributions
imm = MPC.HexIMM()
x = zeros(12)
u = ones(6)
u[1] = 0
y = ones(12)*1e-3
MPC.update!(imm, x, u, y)

(; weights, modes, u_noms, T, obs_dist) = imm
mode = 6
sys = modes[mode]
u_nom = u_noms[mode]
δu = u - u_nom
xp = MPC.dstep(sys, x, δu)
xp
_pdf = pdf(obs_dist(xp), y) # corrector
weights[mode] *= _pdf
