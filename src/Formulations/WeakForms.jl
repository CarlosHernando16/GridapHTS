# ──────────────────────────────────────────────
# Weak Form Building Blocks
# ──────────────────────────────────────────────
#
# Common bilinear forms, linear forms, and residuals used across
# the A-formulation and T-A formulation. Follows Gridap conventions
# for measure-based integration.
# ──────────────────────────────────────────────

"""
    a_bilinear_form(μ_inv, dΩ)

Construct the curl-curl bilinear form for the A-formulation:

    a(A, v) = ∫ (1/μ) (∇×A)⋅(∇×v) dΩ

# Arguments
- `μ_inv`: Inverse permeability 1/μ [m/H] (scalar or function)
- `dΩ`: Domain measure

# Returns
- Bilinear form function `(A, v) -> ...`
"""
function a_bilinear_form(μ_inv, dΩ)
    return (A, v) -> ∫(μ_inv * (∇ × A) ⋅ (∇ × v))dΩ
end

"""
    a_linear_form(f, dΩ)

Construct the source linear form for the A-formulation:

    l(v) = ∫ f⋅v dΩ

where `f` is the source current density J_s.

# Arguments
- `f`: Source function (current density or forcing term)
- `dΩ`: Domain measure

# Returns
- Linear form function `v -> ...`
"""
function a_linear_form(f, dΩ)
    return v -> ∫(f ⋅ v)dΩ
end

"""
    grad_grad_bilinear(dΩ)

Construct the grad-grad bilinear form (Poisson-type):

    a(u, v) = ∫ ∇u⋅∇v dΩ

Useful for validation and scalar potential problems.
"""
function grad_grad_bilinear(dΩ)
    return (u, v) -> ∫(∇(u) ⋅ ∇(v))dΩ
end

"""
    mass_bilinear(dΩ)

Construct the L² mass bilinear form:

    m(u, v) = ∫ u⋅v dΩ
"""
function mass_bilinear(dΩ)
    return (u, v) -> ∫(u ⋅ v)dΩ
end

"""
    ta_residual(mat::AbstractHTSMaterial, μ_inv, dΩ, dΩ_sc)

Construct the residual for the coupled T-A formulation.

The T-A residual consists of two coupled equations:
1. **A-equation** (whole domain): (1/μ)(∇×A)⋅(∇×v_A) = 0
2. **T-equation** (superconductor only): E(J)⋅∇v_T = 0

where J = -∇T (in the thin-tape approximation) and E = ρ(|J|)J.

# Arguments
- `mat::AbstractHTSMaterial`: HTS material model
- `μ_inv`: Inverse permeability
- `dΩ`: Whole-domain measure
- `dΩ_sc`: Superconductor subdomain measure

# Returns
- Residual function `((A, T), (v_A, v_T)) -> ...`
"""
function ta_residual(mat::PowerLawMaterial, μ_inv, dΩ, dΩ_sc)
    return ((A, T), (v_A, v_T)) -> begin
        # A-equation: curl-curl over full domain
        res_A = ∫(μ_inv * (∇ × A) ⋅ (∇ × v_A))dΩ

        # T-equation: nonlinear resistivity in superconductor
        # J = -∇T, E = ρ(|J|) * J
        J = -∇(T)
        J_norm = norm_safe(J)
        rho = resistivity(mat, J_norm)
        res_T = ∫(rho * J ⋅ ∇(v_T))dΩ_sc

        res_A + res_T
    end
end

"""
    ta_jacobian(mat::PowerLawMaterial, μ_inv, dΩ, dΩ_sc)

Construct the Jacobian for the T-A formulation (for Newton-Raphson).

# Arguments
Same as [`ta_residual`](@ref).

# Returns
- Jacobian function `((A, T), (dA, dT), (v_A, v_T)) -> ...`
"""
function ta_jacobian(mat::PowerLawMaterial, μ_inv, dΩ, dΩ_sc)
    return ((A, T), (dA, dT), (v_A, v_T)) -> begin
        # A-equation Jacobian: linear curl-curl (w.r.t. dA)
        jac_AA = ∫(μ_inv * (∇ × dA) ⋅ (∇ × v_A))dΩ

        # T-equation Jacobian: nonlinear part (w.r.t. dT)
        J = -∇(T)
        J_norm = norm_safe(J)
        rho = resistivity(mat, J_norm)
        drho = resistivity_derivative(mat, J_norm)

        # ρ(|J|) * ∇(dT)⋅∇(v_T) + dρ/d|J| * (J⋅∇(dT)/|J|) * J⋅∇(v_T)
        dJ = -∇(dT)
        jac_TT = ∫(
            rho * dJ ⋅ ∇(v_T) +
            drho * ((J ⋅ dJ) / J_norm) * (J ⋅ ∇(v_T))
        )dΩ_sc

        jac_AA + jac_TT
    end
end
