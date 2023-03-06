const m = 2.4 # kg
const g = 9.81
const Jx = 5.126E-3 # kgm^2
const Jy = 5.126E-3 # kgm^2
const Jz = 1.3E-2 # kgm^2
const b = 2.98E-5 # N/rad^2
const l = 0.5 # M
const d = 1.140E-7 # N*m/rad^2

# Define State Space Model
const nn = 12 # num state variables
const nmFM = 4 # num of virtual
const nm = 6 # num actuator inputs
const np = 12 # num outputs

const STATE_LABELS = [
    "x";;
    "y";;
    "z";;
    "dx";;
    "dy";;
    "dz";;
    "ϕ";;
    "θ";;
    "ψ";;
    "dϕ";;
    "dθ";;
    "dψ"
]

const TRANSLATIONAL_STATES = [1,2,3,7,8,9]
const ROTATIONAL_STATES = [4,5,6,10,11,12]
