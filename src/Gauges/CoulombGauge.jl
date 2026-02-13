# ──────────────────────────────────────────────
# Coulomb Gauge: div(A) = 0
# ──────────────────────────────────────────────
#
# The Coulomb gauge ensures uniqueness of the magnetic vector potential A
# by imposing ∇⋅A = 0. Implemented via a Lagrange multiplier approach
# that creates a saddle-point system.
# ──────────────────────────────────────────────

"""
    apply_coulomb_gauge(V_A, U_A, model, Ω, dΩ, order;
        dirichlet_tags=String[]) -> NamedTuple

Augment the A-formulation FE spaces with a Lagrange multiplier λ
to enforce the Coulomb gauge condition ∇⋅A = 0.

The resulting saddle-point system is:

    [K   G] [A]   [f]
    [G'  0] [λ] = [0]

where G corresponds to the constraint ∫ (∇⋅A) λ dΩ.

# Arguments
- `V_A`: Test FE space for A (H(curl) or H¹)
- `U_A`: Trial FE space for A
- `model`: Discrete model
- `Ω`: Triangulation
- `dΩ`: Measure
- `order::Int`: Polynomial order for Lagrange multiplier space
- `dirichlet_tags`: Dirichlet boundary tags for λ (default: none)

# Returns
Named tuple with:
- `Y`: MultiField test space [V_A, V_λ]
- `X`: MultiField trial space [U_A, U_λ]
- `constraint_form`: Bilinear form for the gauge constraint

# Theory
The Coulomb gauge div(A) = 0 is imposed weakly via:

    ∫ (∇⋅A) λ dΩ = 0  ∀λ

This creates a saddle-point structure similar to the pressure in
Stokes flow.

# References
- Creusé et al. (2018), Math. Methods Appl. Sci., 41(16), 6585-6612.
"""
function apply_coulomb_gauge(V_A, U_A, model, Ω, dΩ, order;
    dirichlet_tags=String[]
)
    # Lagrange multiplier space (H¹, scalar)
    reffe_λ = ReferenceFE(lagrangian, Float64, order)

    if isempty(dirichlet_tags)
        V_λ = TestFESpace(model, reffe_λ; conformity=:H1)
    else
        V_λ = TestFESpace(model, reffe_λ;
            conformity=:H1, dirichlet_tags=dirichlet_tags)
    end
    U_λ = TrialFESpace(V_λ, 0.0)

    # MultiField spaces
    Y = MultiFieldFESpace([V_A, V_λ])
    X = MultiFieldFESpace([U_A, U_λ])

    # Gauge constraint bilinear form: ∫ (∇⋅A) λ dΩ
    constraint_form = (A, λ) -> ∫(divergence(A) * λ)dΩ

    return (
        Y = Y,
        X = X,
        V_λ = V_λ,
        U_λ = U_λ,
        constraint_form = constraint_form,
    )
end
