# sin(30) = 0.5
# cos(30) = √3 / 2

const s_30 = sind(30)
const c_30 = cosd(30)

const HOVER_LHS = @SMatrix [
    -s_30 -1 -s_30 s_30 1 s_30 # roll
    c_30 0 c_30 -c_30 0 -c_30 # pitch
    1 -1 1 -1 1 -1 # yaw
    1 1 1 1 1 1
]

const HOVER_RHS = @SVector [0,0,0,m*g]

function hover_control(failed=0)
    if !iszero(failed)
        not_failed = filter(!=(failed), 1:6)
        A = HOVER_LHS[:,not_failed]
        u_not_failed = pinv(A)*HOVER_RHS
        u_eq = zeros(6)
        u_eq[not_failed] .= u_not_failed
        return u_eq
    else
        return pinv(HOVER_LHS)*HOVER_RHS
    end
end
