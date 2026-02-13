# ──────────────────────────────────────────────
# Neumann (Natural) Boundary Conditions
# ──────────────────────────────────────────────

"""
    apply_neumann_bc(g, dΓ)

Construct the Neumann boundary contribution to the linear form:

    l_N(v) = ∫ g⋅v dΓ

where `g` is the prescribed boundary flux/traction and `dΓ` is the
boundary measure.

# Arguments
- `g`: Neumann data (function or constant)
- `dΓ`: Boundary measure (from `Measure(Γ, degree)`)

# Returns
- Linear form contribution `v -> ∫ g⋅v dΓ`

# Example
```julia
# Applied magnetic field on external boundary
Γ = BoundaryTriangulation(model, tags=["outer"])
dΓ = Measure(Γ, 2*order)
g(x) = VectorValue(0.0, B_applied)
l_neumann = apply_neumann_bc(g, dΓ)
```
"""
function apply_neumann_bc(g, dΓ)
    return v -> ∫(g ⋅ v)dΓ
end

"""
    apply_neumann_bc(g::Real, dΓ)

Scalar Neumann condition (e.g., for 2D A-formulation).
"""
function apply_neumann_bc(g::Real, dΓ)
    return v -> ∫(g * v)dΓ
end
