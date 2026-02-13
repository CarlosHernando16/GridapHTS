# GridapHTS.jl

*Finite element solver for High-Temperature Superconductor electromagnetics*

## Overview

GridapHTS implements the **T-A formulation** (coupled magnetic vector potential **A** and electric scalar potential **T**) for modeling HTS screening currents, nonlinear power-law J-E behavior, and quench dynamics.

Built on the [Gridap.jl](https://github.com/gridap/Gridap.jl) ecosystem, following [GridapMHD.jl](https://github.com/gridapapps/GridapMHD.jl) architectural patterns.

## Key Features

- **A-formulation**: Magnetic vector potential for magnetostatic problems
- **T-A formulation**: Coupled A and T for HTS current distribution
- **Nonlinear power-law**: ``\mathbf{E} = E_c (J/J_c)^n \hat{\mathbf{J}}`` constitutive relation
- **Field-dependent ``J_c``**: Kim model and custom ``J_c(B,T)`` dependencies
- **Newton-Raphson with n-continuation**: Robust nonlinear solver strategy
- **Coulomb gauge**: ``\nabla \cdot \mathbf{A} = 0`` via Lagrange multiplier
- **DrWatson workflow**: Reproducible simulations with parameter sweeps

## Quick Start

```julia
using DrWatson
@quickactivate "GridapHTS"

using GridapHTS

params = Dict(
    :formulation => :A,
    :mesh => Dict(:type => :cartesian, :domain => (0,1,0,1), :partition => (20,20)),
    :fe_order => 1,
)

xh, fullparams, info = GridapHTS.main(params)
```

## Development Phases

| Phase | Description | Status |
|:------|:------------|:-------|
| 1 | Foundation & 2D validation (A-formulation) | In Progress |
| 2 | Nonlinear T-A with IESL benchmark validation | Planned |
| 3 | 3D geometries + MPI parallelization | Planned |
| 4 | Transient electromagnetics & quench modeling | Planned |

## Contents

```@contents
Pages = ["formulations.md", "benchmarks.md", "api.md"]
Depth = 2
```
