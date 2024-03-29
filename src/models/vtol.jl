const VTOL_Ts = 0.1 # s

VTOL_models() = [fault_free_vtol(), actuator_fault_vtol(), system_fault_vtol()]

function fault_free_vtol()
    F = [
        0.9964 0.0026 -0.0004 -0.0460;
        0.0045 0.9037 -0.0188 -0.3834;
        0.0098 0.0339 0.9383 0.1302;
        0.0005 0.0017 0.0968 1.0067;
    ]

    G = [
        0.0445 0.0167;
        0.3407 -0.7249;
        -0.5278 0.4214;
        -0.0268 0.0215
    ]

    H = [
        1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 1.0 1.0 1.0
    ]

    M = default_d(G)
    return DTLinearModel(ss(F, G, H, M, VTOL_Ts))
end

function actuator_fault_vtol()
    F = [
        0.9964 0.0026 -0.0004 -0.0460;
        0.0045 0.9037 -0.0188 -0.3834;
        0.0098 0.0339 0.9383 0.1302;
        0.0005 0.0017 0.0968 1.0067;
    ]

    G = [
         0.0045  0.0167;
         0.0000 -0.3624;
        -0.1319  0.1053;
        -0.0268  0.0215
    ]

    H = [
        1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 1.0 1.0 1.0
    ]

    M = default_d(G)
    return DTLinearModel(ss(F, G, H, M, VTOL_Ts))
end

function system_fault_vtol()
    F = [
        0.9964 0.0026 -0.0004 -0.0460;
        0.0045 0.0000 -0.0188 -0.3834;
        0.0098 0.0339 0.9383 0.1302;
        0.0005 0.0017 0.0968 1.0067;
    ]

    G = [
        0.0445 0.0167;
        0.3407 -0.7249;
        -0.5278 0.4214;
        -0.0268 0.0215
    ]

    H = [
        1.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 1.0 1.0 1.0
    ]

    M = default_d(G)
    return DTLinearModel(ss(F, G, H, M, VTOL_Ts))
end
