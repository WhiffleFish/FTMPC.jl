function SecondOrderIntegrator(; cx=1.0, cy=1.0)
    A = Float64[
        0 1 0 0;
        0 0 0 0;
        0 0 0 1;
        0 0 0 0
    ]
    B = Float64[
        0 0;
        cx 0;
        0 0;
        0 cy
    ]
    C = default_c(B)
    D = default_d(B)
    return CTLinearModel(ss(A,B,C,D))
end
