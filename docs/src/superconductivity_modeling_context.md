# FEM Superconductivity Modeling Context

This page captures practical modeling guidance for `GridapHTS.jl`, combining:

- FEM and implementation patterns from `arXiv:1910.01412v2` (Gridap user guide),
- Domain knowledge from *Numerical Modeling of Superconducting Applications* (Dutoit, Grilli, Sirois; 2022/2023 OA edition),
- Established HTS electromagnetic modeling practice for power-law and critical-state behavior.

## 1. Physics Scope and Assumptions

For most HTS devices in power and magnet applications, quasimagnetostatic Maxwell equations are used:

- `\nabla \times \mathbf{H} = \mathbf{J}` (displacement current neglected),
- `\nabla \cdot \mathbf{B} = 0`,
- `\nabla \times \mathbf{E} = -\partial_t \mathbf{B}`.

Typical assumptions in this project:

- Non-magnetic media unless specified (`\mu \approx \mu_0`),
- Slowly varying fields (AC/DC ramps, no wave propagation),
- Macroscopic constitutive law `\mathbf{E}(\mathbf{J}, \mathbf{B}, T)` for HTS.

## 2. HTS Constitutive Laws

### Power-law `E-J` model (default)

The standard nonlinear law is:

```math
\mathbf{E} = E_c \left(\frac{|\mathbf{J}|}{J_c(B,T,\theta)}\right)^n \frac{\mathbf{J}}{|\mathbf{J}|}
```

Practical notes:

- `J_c` can be field-, angle-, and temperature-dependent.
- Large `n` approximates the critical-state model.
- Numerical regularization is needed near `|\mathbf{J}| \to 0` (safe norm).

### Critical-state interpretation

The critical-state model remains essential for interpreting penetration fronts, hysteresis, and AC loss trends even when solving the power-law model numerically.

## 3. Formulations Relevant to GridapHTS

### A-formulation

- Strong form: `\nabla \times (\mu^{-1}\nabla \times \mathbf{A}) = \mathbf{J}`.
- In 2D (`A_z` scalar), this becomes a Poisson-type problem.
- In 3D, enforce gauge (typically Coulomb) for uniqueness.

### T-A formulation

Common for coated conductors and thin superconducting regions:

- `\mathbf{A}` equation on full domain,
- `T` equation on superconducting region with `\mathbf{J}` derived from `T`,
- Coupling through source current terms and nonlinear `E-J`.

### Other formulations (when needed)

- `H`-formulation: robust and general but computationally heavier,
- `H-\phi` and reduced formulations: useful for speedups in specific geometries.

## 4. FEM Discretization Guidance

For this codebase and Gridap-style implementation:

- Use conforming `H1` spaces for scalar potentials (`A_z`, `T`, multipliers),
- Use `H(curl)`-conforming spaces for full 3D vector `\mathbf{A}`,
- Keep integration degree high enough for nonlinear terms (`>= 2*order+1`, often `2*order+2` for power-law terms),
- Apply physically consistent Dirichlet/Neumann boundaries and transport-current constraints,
- Maintain clear domain tagging for superconductors, air, and stabilizers.

## 5. Nonlinear and Transient Solver Strategy

Recommended workflow:

1. Start from lower `n` (continuation), then ramp to target `n`.
2. Use Newton (or damped Newton) with robust line search.
3. For transients, start with conservative time steps and adapt.
4. Monitor residual, Jacobian conditioning, and energy/loss consistency.
5. Regularize divisions by current magnitude (`norm_safe` pattern).

## 6. Verification and Validation Checklist

Minimum checks for each new formulation or feature:

- Recovery of known magnetostatic limits in linear cases,
- Current penetration and hysteresis trends versus analytical/CSM expectations,
- AC loss comparison versus benchmark geometries (single tape, stacks, coils),
- Mesh and time-step convergence study,
- Sensitivity to `J_c(B,T,\theta)` parameterization,
- Cross-check with published reference curves where available.

## 7. Multiphysics Context (Roadmap)

The 2022 superconducting modeling reference emphasizes integrated simulation:

- Electromagnetics,
- Thermal stability and quench,
- Thermo-hydraulics,
- Structural/mechanical effects.

For `GridapHTS.jl`, this motivates clean interfaces so EM solvers can later couple to thermal/mechanical modules without redesigning core unknowns and material APIs.

## 8. Implementation Patterns from Gridap User Guide

The Gridap paper (`1910.01412v2`) is not HTS-specific but is highly relevant to project architecture:

- Express weak forms close to mathematical notation,
- Build FE spaces, triangulations, quadrature, and terms explicitly,
- Separate operator construction from solver strategy,
- Preserve composability for multifield and nonlinear operators.

These patterns should guide all new formulations in this repository.

## References

- F. Verdugo, S. Badia, *A User-Guide to Gridap* (`arXiv:1910.01412v2`).
- B. Dutoit, F. Grilli, F. Sirois (eds.), *Numerical Modeling of Superconducting Applications* (World Scientific, OA edition, 2022/2023).
