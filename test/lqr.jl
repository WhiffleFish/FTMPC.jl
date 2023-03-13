x0 = zeros(12)
x_ref = zeros(12)
x_ref[3] = -5
Q = I; R = I
t = 0:0.1:100

for mode in 0:6
    model = LinearHexModel(mode)
    L = lqr(model.ss, Q, R)
    u(x,t) = -L*(x - x_ref)
    y, t, x, uout = lsim(model.ss,u,t,x0=x0)
    @test xout ≈ x_ref atol=1e-3
    @test uout ≈ zeros(6) atol=1e-3
end
