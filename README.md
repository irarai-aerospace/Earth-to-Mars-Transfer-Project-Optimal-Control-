# Earth-to-Mars-Transfer-Project-Optimal-Control-
# Earth-to-Mars Transfer Project — Optimal Control

A MATLAB codebase covering orbital mechanics, trajectory optimization, stochastic propagation, and convex-programming-based landing for interplanetary and Earth-orbit missions. The scripts span classical and modern orbital element representations, indirect optimal control (Pontryagin Minimum Principle), convex relaxation, and linear/Monte Carlo uncertainty propagation.

---

## File Summaries

### Coordinate & Orbital Element Conversion Utilities

---

#### `keplerian_to_cartesian.m`
**Concept:** Converts classical Keplerian orbital elements to an inertial Cartesian state vector by solving Kepler's equation for the eccentric anomaly and applying a 3-2-3 Euler rotation sequence.  
**Inputs:** Semi-major axis `a`, eccentricity `e`, mean anomaly `M`, inclination `i`, RAAN `Ω` (sig), argument of periapsis `ω` (om), gravitational parameter `μ`.  
**Outputs:** Inertial position vector `R_inertial` (3×1) and velocity vector `V_inertial` (3×1).

---

#### `equinoctial_to_cartesian.m`
**Concept:** Converts modified equinoctial orbital elements (MEE) — which avoid singularities near zero eccentricity and zero/180° inclination — to a Cartesian state. Internally re-derives classical elements then delegates to `keplerian_to_cartesian`.  
**Inputs:** MEE components `p, f, g, h, k, L` (semi-latus rectum, equinoctial eccentricity and inclination projections, and true longitude), gravitational parameter `μ`.  
**Outputs:** Inertial position `R_inertial` (3×1) and velocity `V_inertial` (3×1).

---

#### `milankovitch_to_cartesian.m`
**Concept:** Converts Milankovitch orbital elements — the angular momentum vector **h**, eccentricity vector **e**, and true longitude `L` — to a Cartesian state. This non-singular representation is well suited for nearly circular or equatorial orbits.  
**Inputs:** Angular momentum vector `h` (3×1), eccentricity vector `e` (3×1), true longitude `L` (rad), gravitational parameter `μ`.  
**Outputs:** Inertial position `R_inertial` (3×1) and velocity `V_inertial` (3×1).

---

### Orbit Simulation Scripts

---

#### `keplerian_orbit_simulation.m`
**Concept:** Propagates a Keplerian orbit using the Gauss variational equations in classical element form (`[a, e, i, Ω, ω, M]`) with an adaptive RK4 integrator over 15 orbital periods. In the unperturbed two-body case only mean motion (`Ṁ = n`) is nonzero, giving a steady circular/elliptic orbit.  
**Inputs (hardcoded):** `a₀ = 9000 km`, `e = 0.2`, `i = 50°`, `Ω = 10°`, `ω = 20°`, `μ = 3.986×10¹⁴ m³/s²`.  
**Outputs:** 3D trajectory plot of the orbit over 15 revolutions; also prints the array dimensions of the integrated state history.

---

#### `equinoctial_orbit_simulation.m`
**Concept:** Same physical orbit as `keplerian_orbit_simulation.m` but expressed in modified equinoctial elements. Propagates the Gauss variational equations for MEE under pure two-body motion, then converts each time step back to Cartesian for a 3D trajectory plot.  
**Inputs (hardcoded):** Same orbital parameters as above, converted to MEE at initialization. Adaptive RK4 with `RelTol = AbsTol = 1×10⁻⁶`.  
**Outputs:** 3D orbit trajectory plot (axis-equal, grid on).

---

### Discretization & State Transition Matrix Tools

---

#### `discretization.m`
**Concept:** Derives the zero-order-hold (ZOH) discrete-time matrices for the Mars-landing dynamical system using a matrix exponential for `Aₖ` and numerical (trapezoidal) integration for `Bₖ` and `cₖ`. This is the first step needed before any discrete convex optimization.  
**Inputs (hardcoded):** Time step `dt = 1 s`, thrust-specific-impulse coefficient `α = 0.5086 s/km`, Mars gravity vector `g = [−3.71×10⁻³; 0; 0] km/s²`. Continuous-time matrices `A`, `B`, `c` for the 7-state landing model (`[r; v; ln(m)]`).  
**Outputs:** Printed discrete matrices `Aₖ` (7×7), `Bₖ` (7×4), `cₖ` (7×1) to the MATLAB console.

---

#### `stm_p2.m`
**Concept:** Integrates the 6-DOF two-body + J2 perturbation equations of motion together with the State Transition Matrix (STM) Φ(t, t₀), whose dynamics are `Φ̇ = A(x)Φ`, where `A = ∂f/∂x` is computed symbolically via MATLAB's Symbolic Toolbox. This enables linear sensitivity analysis around a reference trajectory.  
**Inputs (hardcoded):** `a = 9000 km`, `e = 0.2`, `i = 50°`, `Ω = 10°`, `ω = 20°`, `M₀ = 20°`, `μ = 3.986×10¹⁴ m³/s²`, `J₂ = 1.08263×10⁻³`, `Rₑ = 6378 km`. Integration over `tf = 15` orbital periods.  
**Outputs:** Final STM `Φ(tf, t₀)` (6×6 matrix) printed to console; intermediate trajectory and STM history are available in memory.

---

### Optimal Control — Indirect Shooting Methods (Earth-to-Mars Transfer)

All three shooting scripts use a planar heliocentric model with nondimensional units (μ = 1, 1 AU length scale). Earth starts at `a = 1 AU`, `ν₀ = 0`; Mars starts at `a = 1.524 AU`, `ν₀ = π`. Transfer time is `tf = 8` (nondimensional) unless noted.

---

#### `equinoctial_transfer_min_energy.m`  &  `shooting_method_orbit_transfer.m`
**Concept (identical):** Solves the Earth-to-Mars **minimum-energy** (minimum ∫‖u‖² dt) fixed-time transfer using single shooting. The Pontryagin Minimum Principle yields `u*(t) = −½ Bᵀλ(t)`. The unknown initial costate `λ₀ ∈ ℝ⁴` is found by solving the boundary-value problem `x(tf; λ₀) = x_Mars(tf)` with `fsolve`. Twenty random initial guesses are tried to locate solutions.  
**Inputs:** Earth/Mars initial conditions generated internally; initial costate guesses drawn from `U(−10, 10)⁴`.  
**Outputs:** Converged `λ₀`, spacecraft trajectory, optimal control history, costate history, and Hamiltonian evolution `H(t) − H(0)` (should be ≈ 0) — all plotted. `shooting_method_orbit_transfer.m` also stores per-run convergence history.

---

#### `min_fuel_optimal_trajectory.m`
**Concept:** Solves the Earth-to-Mars **minimum-fuel** (minimum ∫γ dt) fixed-time transfer. The optimal control structure is bang-bang with `‖u‖ ≤ u_max = 0.1`. The hard switching function is smoothed via a `tanh` homotopy: `γ = ½ u_max (1 + tanh(S/ρ))`, where `S = ‖Bᵀλ‖ − 1`. A **ρ-continuation** sequence `ρ ∈ {1, 0.1, 0.01, 0.001}` progressively sharpens toward the true bang-bang law.  
**Inputs:** Fixed initial costate guess `λ₀ = [0.8; 0.1; 0.2; 1.1]`; continuation parameters as above.  
**Outputs:** Optimal trajectory, control history, throttle `Γ(t)`, switching function `S(t)`, costate history, Hamiltonian conservation check. Also **saves** `ps4_reference_traj.mat` (fields: `t_ref`, `x_ref`, `u_ref`, `tf`, `mu`, `u_max`) for use in PS10.

---

#### `min_fuel_random_intitial_costate.m`
**Concept:** A robustness study that repeats the minimum-fuel shooting with N = 20 random initial costate guesses (`U(−10, 10)⁴`), each with the full ρ-continuation sequence, to assess how sensitive convergence is to initialization (expected to be much harder than the minimum-energy or minimum-time cases).  
**Inputs:** 20 random `λ₀` draws; same Earth/Mars conditions and continuation parameters as `min_fuel_optimal_trajectory.m`.  
**Outputs:** Per-run convergence summary table (exitflag, total iterations, terminal residual norm, failure location in ρ); stem plots of convergence outcomes and residual norms; one representative successful trajectory with throttle and switching function plots.

---

#### `min_time_optimal_transfer_uconst.m`
**Concept:** Solves the Earth-to-Mars **minimum-time** transfer with free final time `tf`. The unknown vector is `y = [λ₀; tf] ∈ ℝ⁵`, with five conditions: state matching `x(tf) = x_Mars(tf)` (4 eqs) and the transversality condition `H(tf) = 0` (1 eq). The optimal control saturates at `u*(t) = −u_max · Bᵀλ / ‖Bᵀλ‖`.  
**Inputs:** Initial costate/time guess `[5; 2; 2; 7; 3]`; `u_max = 0.1`.  
**Outputs:** Optimal `λ₀`, `tf`, spacecraft trajectory vs. Mars trajectory, control history, costate history, and Hamiltonian evolution `H(t) − H(0)` plotted.

---

### Stochastic Propagation & Uncertainty Quantification

---

#### `PS9.m`
**Concept:** Two-part Monte Carlo study of stochastic processes. Part 1 generates M = 20 sample paths of Brownian increments `Δw₁,ₖ` and `Δw₂,ₖ` (i.i.d. Gaussian, `σ = √dt`) and plots them with theoretical 3σ bounds. Part 2 propagates a 2D linear stochastic system `xₖ₊₁ = xₖ + G Δwₖ` (diffusion matrix `G = diag(σ₁, σ₂)`) to build Monte Carlo state histories with expanding 3σ bounds.  
**Inputs (hardcoded):** `dt = 10⁻³`, `T = 1`, M = 20 samples, `σ₁ = 2`, `σ₂ = 3`.  
**Outputs:** Plots of Brownian increment paths and 2D state paths with 3σ envelopes; console printout of empirical vs. theoretical standard deviations.

---

#### `ps10_b.m`
**Concept (Problem Set 10, Part b):** Applies a Brownian velocity disturbance to the optimal minimum-fuel Earth-to-Mars trajectory (loaded from `ps4_reference_traj.mat`) in a Monte Carlo fashion. Each of M = 20 sample trajectories starts from a perturbed initial state (drawn from `P₀ = diag(σ_pos², σ_pos², σ_vel², σ_vel²)`) and is integrated with additive Brownian disturbances `G dw` at each step using a simple Euler-Maruyama scheme.  
**Inputs:** `ps4_reference_traj.mat` (reference trajectory and control), `dt = 0.01`, `σ_pos = σ_vel = 10⁻²`, `σ_dist = 10⁻³`, M = 20.  
**Outputs:** Time histories of all four states (rₓ, r_y, vₓ, v_y) for each MC sample vs. the reference; 2D position trajectory plot of MC samples vs. reference.

---

#### `ps10_c.m`
**Concept (Problem Set 10, Part c):** Extends Part b with **linear covariance propagation** using the discrete State Transition Matrix (STM). At each time step the continuous-time Jacobian `A(x̄, ū)` is computed around the reference, and both the discrete transition matrix `Φₖ` and the process-noise covariance `Σₖ = ∫ₜₖᵗₖ₊₁ e^{A(t−τ)} GGᵀ eᵀ^{A(t−τ)} dτ` are obtained via ODE integration. The covariance is propagated as `Pₖ₊₁ = ΦₖPₖΦₖᵀ + Σₖ`. Monte Carlo deviations are compared against the resulting 3σ bounds.  
**Inputs:** Same reference trajectory file; `dt = 0.1` for covariance propagation grid; same noise parameters.  
**Outputs:** Per-state deviation plots showing M = 20 MC sample deviations overlaid with ±3σ linear covariance bounds.

---

### Convex Optimization

---

#### `cvx.m`
**Concept:** A problem-set exercise file demonstrating how to formulate and solve several convex optimization problems with the CVX disciplined convex programming toolbox. Problems include minimum-‖x‖₂ under a linear constraint, minimum-sum-of-squares under the same constraint, optimal vector for a matrix residual, and symmetric positive-semidefinite / negative-semidefinite / box-constrained matrix approximation problems.  
**Inputs (hardcoded):** Matrices `A` (5×10), `B` (5×5), vector `y` (5×1).  
**Outputs:** Optimal variables printed/stored by CVX; no plots.

---

#### `convex_ars_landing.m`
**Concept (PS6 — Powered Descent Guidance):** Solves a Mars powered-descent (rocket-powered landing) problem using **lossless convexification**. The nonconvex minimum-fuel problem with thrust-magnitude bounds is lifted to a SOCP (second-order cone program) by introducing a slack variable `σ` for thrust magnitude. Constraints include: glide-slope cone `‖[r₂; r₃]‖ ≤ cot(γ_min) r₁`, non-negative altitude `r₁ ≥ 0`, lateral displacement limit `‖[r₂; r₃]‖ ≤ ρ_max`, log-mass lower bound (mass depletion by thrust), and first-order Taylor-linearized upper thrust bound to preserve convexity. Solved with CVX.  
**Inputs (hardcoded):** `N = 78` intervals, `dt = 1 s`, initial position `r₀ = [1.5; 0; 2.0] km`, initial velocity `v₀ = [−0.075; 0.03; 0.1] km/s`, initial mass `m₀ = 2000 kg`, `T_min = 4.97 kN`, `T_max = 13.26 kN`, glide-slope `γ_min = 3°`.  
**Outputs:** Optimal fuel consumption and final mass printed to console; plots of 3D trajectory, position/velocity/mass time histories, lossless convexification check (`σ` vs. `‖a‖`), and log-mass state with lower bound.

---

## Dependencies & Requirements

- **MATLAB** (R2019b or later recommended)
- **CVX** toolbox — required by `cvx.m` and `convex_ars_landing.m`
- **Symbolic Math Toolbox** — required by `stm_p2.m`
- Helper functions expected in the path: `Kepler_equation`, `rk4_adaptive`, `rk4singlestep`, `cartesian_to_equinoctial`

## Data Files

| File | Produced by | Consumed by |
|------|-------------|-------------|
| `ps4_reference_traj.mat` | `min_fuel_optimal_trajectory.m` | `ps10_b.m`, `ps10_c.m` |
