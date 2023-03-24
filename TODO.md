# TODO

- [x] Hex Linear
  - [x] Fully Observable
    - [x] LQR
    - [x] OSQP LQR
    - [x] JuMP
- [x] Support multiple CBFs/constraints
- [ ] Simulation
  - [x] Parameter Updater
  - [x] Variable mode weighting
  - [ ] Better warm-starting
    - Use SHIFT of last solution rather than last solution directly
  - [x] Variable Consensus Horizon
- [ ] IMM integration
- [ ] Feasibility LP (No QP)
- [ ] Change `LinearHexModel` to discrete time (we're always converting anyways)
