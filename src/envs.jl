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

function ElevatorShaft(;w=1,l=1,h=3,γg=1e-1, γs=1e-1)
    return [
        LinearConstraint(-basis(12, 3)*1, h, γg),
        LinearConstraint(basis(12, 2)*1, l, γs),
        LinearConstraint(-basis(12, 2)*1, l, γs),
        LinearConstraint(basis(12, 1)*1, w, γs),
        LinearConstraint(-basis(12, 1)*1, w, γs)
    ]
end


function RendezvousBarrier(;w=1,l_lower=2,l_upper=5,h=3,γg=1e-1, γs=1e-1)
    return [
        LinearConstraint(basis(6, 3)*1, h, γg),
        LinearConstraint(-basis(6, 3)*1, h, γg),
        LinearConstraint(-basis(6, 2)*1, l_upper, γs),
        LinearConstraint(basis(6, 2)*1, -l_lower, γs),
        LinearConstraint(basis(6, 1)*1, w, γs),
        LinearConstraint(-basis(6, 1)*1, w, γs)
    ]
end

# LinearConstraint(basis(6, 2)*1, -2, γs) -> y>=2
# LinearConstraint(basis(6, 2)*1, 5, γs) -> y<=5
