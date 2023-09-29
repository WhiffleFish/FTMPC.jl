function Box(;w=1,l=1,h=1,γ=1e-1)
    return [
        LinearConstraint(basis(12, 3)*1, h, γ),
        LinearConstraint(-basis(12, 3)*1, h, γ),
        LinearConstraint(basis(12, 2)*1, l, γ),
        LinearConstraint(-basis(12, 2)*1, l, γ),
        LinearConstraint(basis(12, 1)*1, w, γ),
        LinearConstraint(-basis(12, 1)*1, w, γ)
    ]
end

function ElevatorShaft(;w=1,l=1,h=3,γ=1e-1)
    return [
        LinearConstraint(-basis(12, 3)*1, h, γ),
        LinearConstraint(basis(12, 2)*1, l, γ),
        LinearConstraint(-basis(12, 2)*1, l, γ),
        LinearConstraint(basis(12, 1)*1, w, γ),
        LinearConstraint(-basis(12, 1)*1, w, γ)
    ]
end

function VTOLCeiling(sys::BatchDynamics; h=1, γ=1e-1)
    return [
        LinearConstraint(basis(inner_statedim(sys), 2)*1, h, γ),
        LinearConstraint(-basis(inner_statedim(sys), 2)*1, h, γ)
    ]
end