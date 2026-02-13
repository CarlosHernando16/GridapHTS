# Getting Started

## Prerequisites

- **Julia >= 1.10**: [julialang.org/downloads](https://julialang.org/downloads/)
- **Git**: Version control
- **Gmsh** (optional): Mesh generation ([gmsh.info](https://gmsh.info/))
- **ParaView** (optional): Visualization ([paraview.org](https://www.paraview.org/))

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/magnestar/GridapHTS.jl")
```

## First Simulation

```julia
using DrWatson
@quickactivate "GridapHTS"

using GridapHTS

# Solve a simple magnetostatic problem
params = Dict(
    :formulation => :A,
    :fe_order => 1,
    :mesh => Dict(
        :type => :cartesian,
        :domain => (0, 1, 0, 1),
        :partition => (20, 20),
    ),
    :output => Dict(:vtk => true, :vtk_path => plotsdir("tutorial", "first_sim")),
)

xh, fullparams, info = GridapHTS.main(params)
@info "Done! Elapsed: $(info[:elapsed_time])s"
```

## Viewing Results

Open the generated `.vtu` file in ParaView:

```bash
paraview plots/tutorial/first_sim.vtu
```
