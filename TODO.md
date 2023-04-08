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
- [x] IMM integration
- [ ] Feasibility LP (No QP) - doesn't actually seem to be faster
- [ ] Change `LinearHexModel` to discrete time (we're always converting anyways)
- [x] Early sim stopping for barrier violation / infeasible control
- [x] Fixed failure time - don't want failure in simulation to be up to random chance
