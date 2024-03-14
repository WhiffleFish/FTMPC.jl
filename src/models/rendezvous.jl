function RendezvousModel(n=0.056)
    m = 100.0
    A = Float64[
        0 0 0 1 0 0;
        0 0 0 0 1 0;
        0 0 0 0 0 1;
        3*n^2 0 0 0 2*n 0;
        0 0 0 -2*n 0 0;
        0 0 -n^2 0 0 0;

    ]
    B = Float64[
        0 0 0;
        0 0 0;
        0 0 0;
        1/m 0 0;
        0 1/m 0;
        0 0 1/m
    ]
    C = default_c(B)
    D = default_d(B)
    return CTLinearModel(ss(A,B,C,D))
end
