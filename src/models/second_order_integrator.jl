function DoubleIntegrator(; cx=1.0, cy=1.0, ueff=1.0)
    # [x ẋ y ẏ]
    A = Float64[
        0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0
    ]
    B = Float64[
        0 0;
        0 0;
        cx 0;
        0 ueff*cy
    ]
    C = default_c(B)
    D = default_d(B)
    return ss(A,B,C,D)
end
