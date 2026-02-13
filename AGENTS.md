# GridapHTS AI Agent Instructions

## Project Overview

**GridapHTS** is an open-source Julia finite element solver for High-Temperature Superconductor (HTS) electromagnetics. It implements the **T-A formulation** (coupled magnetic vector potential **A** and electric scalar potential **T**) to model HTS screening currents, nonlinear power-law J-E behavior, and quench dynamics.

Built on **Gridap.jl** (finite elements), following **GridapMHD.jl** architectural patterns, and managed with **DrWatson.jl** (reproducible scientific workflow).

## Architecture

```
src/
├── GridapHTS.jl          # Main module: imports, includes, exports
├── constants.jl          # Physical constants (MU_0, E_C_DEFAULT, etc.)
├── parameters.jl         # add_default_params(), validate_params()
├── utils.jl              # norm_safe(), get_nested(), merge_nested()
├── main.jl               # main(params::Dict) entry point
├── exports.jl            # All public API exports
├── Materials/
│   ├── PowerLaw.jl       # PowerLawMaterial: E = Ec(J/Jc)^n
│   ├── FieldDependence.jl # KimModel: Jc(B) = Jc0/(1+|B|/B0)
│   └── MaterialLibrary.jl # rebco_default(), bscco_default()
├── Formulations/
│   ├── WeakForms.jl      # Bilinear forms, residuals, Jacobians
│   ├── AFormulation.jl   # Magnetostatic A-formulation
│   └── TAFormulation.jl  # Coupled T-A for HTS
├── Solvers/
│   └── NewtonRaphson.jl  # solve_newton(), solve_with_continuation()
├── Meshers/
│   └── GmshIO.jl        # load_gmsh_model(), setup_cartesian_model()
├── BoundaryConditions/
│   ├── Dirichlet.jl      # apply_dirichlet_bc()
│   └── Neumann.jl        # apply_neumann_bc()
├── Gauges/
│   └── CoulombGauge.jl   # apply_coulomb_gauge() via Lagrange multiplier
└── Applications/
    └── IESLTape.jl       # setup_iesl_benchmark()
```

## Key Design Patterns

### 1. Parameter-Driven Architecture
Everything flows through a `Dict{Symbol,Any}`:
```julia
params = Dict(
    :formulation => :A,      # or :TA
    :fe_order => 1,
    :mesh => Dict(:type => :cartesian, :domain => (0,1,0,1), :partition => (20,20)),
    :material => Dict(:jc => 1e9, :n_exponent => 25),
    :solver => Dict(:type => :newton, :max_iter => 50, :rtol => 1e-8),
    :output => Dict(:vtk => true),
)
xh, fullparams, info = GridapHTS.main(params)
```

### 2. Val-Based Dispatch (from GridapMHD)
Formulation selection uses `Val` types:
```julia
_solve_formulation(::Val{:A}, params)   # A-formulation
_solve_formulation(::Val{:TA}, params)  # T-A formulation
```

### 3. DrWatson Conventions
- Scripts start with `using DrWatson; @quickactivate "GridapHTS"`
- Paths: `datadir()`, `plotsdir()`, `scriptsdir()`
- Results: `@tagsave`, `produce_or_load`

### 4. Material Abstraction
`AbstractHTSMaterial` -> `PowerLawMaterial`, `FieldDependentMaterial`
Each implements `resistivity(mat, J_norm)` and `electric_field(mat, J)`.

## Common Development Tasks

### Adding a New Material Model
1. Create `src/Materials/YourModel.jl`
2. Define struct subtyping `AbstractHTSMaterial`
3. Implement `resistivity()` and `electric_field()`
4. Add dispatch case in `material_from_params()` (`MaterialLibrary.jl`)
5. `include()` in `src/GridapHTS.jl`
6. Export in `src/exports.jl`
7. Add tests in `test/test_materials.jl`

### Adding a New Formulation
1. Create `src/Formulations/YourFormulation.jl`
2. Implement `setup_X_formulation(params)` returning named tuple
3. Implement `solve_X_formulation(params)` convenience function
4. Add `_solve_formulation(::Val{:X}, params)` dispatch in `main.jl`
5. Add VTK output case in `_write_vtk()` in `main.jl`
6. `include()` in `src/GridapHTS.jl`, export in `exports.jl`
7. Add tests in `test/test_formulations.jl`

### Adding a New Solver
1. Create `src/Solvers/YourSolver.jl`
2. Implement solver function accepting FE operator + solver params
3. `include()` in `src/GridapHTS.jl`, export in `exports.jl`
4. Add tests in `test/test_solvers.jl`

### Adding a New Benchmark Application
1. Create `src/Applications/YourBenchmark.jl`
2. Implement `setup_X_benchmark(; kwargs...)` returning params Dict
3. `include()` in `src/GridapHTS.jl`, export in `exports.jl`
4. Create script in `scripts/phaseN/run_X.jl`

## Physics Reference

### Power-Law E-J Relation
```
E = E_c * (|J| / J_c)^n * (J / |J|)
rho(|J|) = (E_c / J_c) * (|J| / J_c)^(n-1)
```
- E_c = 1e-4 V/m (critical electric field criterion)
- J_c ~ 1e9 to 3e10 A/m² (critical current density)
- n = 15-50 (power-law exponent; higher = sharper transition)

### A-Formulation
```
curl(1/mu * curl(A)) = J_s
```
In 2D: A is scalar (out-of-plane), reduces to scalar Poisson.

### T-A Formulation
```
A-equation: curl(1/mu * curl(A)) = J    (full domain)
T-equation: div(E(J)) = 0               (superconductor)
J = -grad(T)
```

### Coulomb Gauge
div(A) = 0 imposed weakly via Lagrange multiplier for 3D uniqueness.

## Testing

Run tests: `julia --project=. -e "using Pkg; Pkg.test()"`

Tests validate:
- Material models (power-law limits, derivatives, field dependence)
- Formulations (Poisson against analytical, A-formulation zero source)
- Parameters (defaults, validation, nested operations)
- Main driver (end-to-end A-formulation)

## References

- [Gridap.jl](https://github.com/gridap/Gridap.jl) - Julia FE framework
- [GridapMHD.jl](https://github.com/gridapapps/GridapMHD.jl) - Architecture reference
- [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/stable/) - Project management
- [htsmodelling.com](https://htsmodelling.com/) - HTS benchmark problems
- Rhyner (1993), Physica C, 212(3-4), 292-300 - Power-law E-J relation
- Kim et al. (1962), Phys. Rev. Lett., 9(7), 306 - Field-dependent J_c
