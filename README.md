# GridapHTS.jl

**Finite element solver for High-Temperature Superconductor electromagnetics**

GridapHTS implements the **T-A formulation** (coupled magnetic vector potential **A** and electric scalar potential **T**) for modeling HTS screening currents, nonlinear power-law J-E behavior, and quench dynamics. Built on the [Gridap.jl](https://github.com/gridap/Gridap.jl) ecosystem.

## Features

- **A-formulation**: Magnetic vector potential for magnetostatic problems
- **T-A formulation**: Coupled A and T for HTS current distribution
- **Nonlinear power-law**: E = E_c (J/J_c)^n constitutive relation
- **Field-dependent J_c**: Kim model and custom J_c(B,T) dependencies
- **Newton-Raphson with n-continuation**: Robust nonlinear solver strategy
- **Coulomb gauge**: div(A) = 0 via Lagrange multiplier
- **Gmsh integration**: Flexible mesh generation via GridapGmsh
- **DrWatson workflow**: Reproducible simulations with parameter sweeps

## Development Phases

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Foundation & 2D validation (A-formulation) | In Progress |
| 2 | Nonlinear T-A with IESL benchmark validation | Planned |
| 3 | 3D geometries + MPI parallelization | Planned |
| 4 | Transient electromagnetics & quench modeling | Planned |

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/magnestar/GridapHTS.jl")
```

Or for development:

```bash
git clone https://github.com/magnestar/GridapHTS.jl.git
cd GridapHTS.jl
julia --project=.
```

```julia
using Pkg
Pkg.instantiate()
```

## Quick Start

```julia
using DrWatson
@quickactivate "GridapHTS"

using GridapHTS

# Define simulation parameters
params = Dict(
    :mesh => Dict(
        :type => :cartesian,
        :domain => (0, 1, 0, 1),
        :partition => (20, 20),
    ),
    :material => Dict(
        :jc => 1e9,
        :n_exponent => 25,
        :ec => 1e-4,
    ),
    :solver => Dict(
        :type => :newton,
        :max_iter => 50,
        :rtol => 1e-8,
    ),
    :formulation => :A,
    :fe_order => 1,
    :vtk_output => true,
)

# Run simulation
xh, fullparams, info = GridapHTS.main(params)
```

## Project Structure

This project follows [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/stable/) conventions:

```
gridapHTS/
├── data/          # Immutable data (meshes, material properties, experimental)
├── docs/          # Documentation (Documenter.jl)
├── meshes/        # Gmsh .geo generation scripts
├── notebooks/     # Jupyter/Pluto notebooks
├── papers/        # Paper drafts & figures
├── plots/         # Generated plots (gitignored)
├── scripts/       # Simulation drivers
├── src/           # Core library source code
└── test/          # Unit & integration tests
```

## References

- [Gridap.jl](https://github.com/gridap/Gridap.jl) - The Julia FE framework
- [GridapMHD.jl](https://github.com/gridapapps/GridapMHD.jl) - Architectural reference
- [DrWatson.jl](https://juliadynamics.github.io/DrWatson.jl/stable/) - Scientific project management
- [htsmodelling.com](https://htsmodelling.com/) - HTS benchmark problems

## License

MIT License. See [LICENSE](LICENSE) for details.
