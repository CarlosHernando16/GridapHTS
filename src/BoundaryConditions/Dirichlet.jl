# ──────────────────────────────────────────────
# Dirichlet (Essential) Boundary Conditions
# ──────────────────────────────────────────────

"""
    apply_dirichlet_bc(V, model, tags, values) -> TrialFESpace

Create a TrialFESpace with Dirichlet boundary conditions applied.

# Arguments
- `V`: Test FE space
- `model`: Discrete model
- `tags`: Boundary tags (String or Vector{String})
- `values`: Boundary values (scalar, vector, or function)

# Returns
- `U::TrialFESpace`: Trial space with Dirichlet BCs

# Example
```julia
# Zero Dirichlet on all boundaries
U = apply_dirichlet_bc(V, model, ["boundary"], 0.0)

# Inhomogeneous Dirichlet
U = apply_dirichlet_bc(V, model, ["inlet"], x -> sin(π*x[2]))
```
"""
function apply_dirichlet_bc(V, model, tags, values)
    return TrialFESpace(V, values)
end

"""
    apply_dirichlet_bc(V, model, tags) -> TrialFESpace

Zero Dirichlet boundary condition (homogeneous).
"""
function apply_dirichlet_bc(V, model, tags)
    return TrialFESpace(V, 0.0)
end

"""
    zero_bc_value(::Val{2}) -> Float64
    zero_bc_value(::Val{3}) -> VectorValue

Return the appropriate zero boundary value for the spatial dimension.
- 2D scalar A: `0.0`
- 3D vector A: `VectorValue(0.0, 0.0, 0.0)`
"""
zero_bc_value(::Val{2}) = 0.0
zero_bc_value(::Val{3}) = VectorValue(0.0, 0.0, 0.0)
