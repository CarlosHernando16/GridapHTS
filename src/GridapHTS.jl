"""
    GridapHTS

Finite element solver for High-Temperature Superconductor (HTS) electromagnetics.

Implements the T-A formulation (coupled magnetic vector potential **A** and
electric scalar potential **T**) for modeling HTS screening currents,
nonlinear power-law J-E behavior, and quench dynamics.

Built on the [Gridap.jl](https://github.com/gridap/Gridap.jl) ecosystem,
following [GridapMHD.jl](https://github.com/gridapapps/GridapMHD.jl) patterns.

# Development Phases
1. Foundation & 2D validation (A-formulation)
2. Nonlinear T-A with IESL benchmark validation
3. 3D geometries + MPI parallelization
4. Transient electromagnetics & quench modeling
"""
module GridapHTS

# ──────────────────────────────────────────────
# Dependencies
# ──────────────────────────────────────────────
using Gridap
using Gridap.FESpaces
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays

using LinearAlgebra
using SparseArrays
using ForwardDiff
using WriteVTK

# ──────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────
include("constants.jl")

# ──────────────────────────────────────────────
# Utilities & Parameters
# ──────────────────────────────────────────────
include("utils.jl")
include("parameters.jl")

# ──────────────────────────────────────────────
# Materials
# ──────────────────────────────────────────────
include("Materials/PowerLaw.jl")
include("Materials/FieldDependence.jl")
include("Materials/MaterialLibrary.jl")

# ──────────────────────────────────────────────
# Mesh Utilities
# ──────────────────────────────────────────────
include("Meshers/GmshIO.jl")

# ──────────────────────────────────────────────
# Boundary Conditions
# ──────────────────────────────────────────────
include("BoundaryConditions/Dirichlet.jl")
include("BoundaryConditions/Neumann.jl")

# ──────────────────────────────────────────────
# Gauge Conditions
# ──────────────────────────────────────────────
include("Gauges/CoulombGauge.jl")

# ──────────────────────────────────────────────
# Formulations (weak forms)
# ──────────────────────────────────────────────
include("Formulations/WeakForms.jl")
include("Formulations/AFormulation.jl")
include("Formulations/TAFormulation.jl")

# ──────────────────────────────────────────────
# Solvers
# ──────────────────────────────────────────────
include("Solvers/NewtonRaphson.jl")

# ──────────────────────────────────────────────
# Main Driver
# ──────────────────────────────────────────────
include("main.jl")

# ──────────────────────────────────────────────
# Applications (benchmark setups)
# ──────────────────────────────────────────────
include("Applications/IESLTape.jl")

# ──────────────────────────────────────────────
# Public API
# ──────────────────────────────────────────────
include("exports.jl")

end # module GridapHTS
