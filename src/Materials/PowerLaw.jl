# ──────────────────────────────────────────────
# Power-Law Constitutive Relation for HTS
# ──────────────────────────────────────────────
#
# The standard E-J power-law for HTS:
#   E = E_c * (|J| / J_c)^n * (J / |J|)
#
# Equivalently, the nonlinear resistivity:
#   ρ(J) = (E_c / J_c) * (|J| / J_c)^(n-1)
# ──────────────────────────────────────────────

"""
    AbstractHTSMaterial

Abstract supertype for all HTS material models.

Subtypes must implement:
- `resistivity(mat, J_norm)` -> scalar resistivity ρ
- `electric_field(mat, J)` -> electric field vector E
"""
abstract type AbstractHTSMaterial end

"""
    PowerLawMaterial <: AbstractHTSMaterial

Standard E-J power-law constitutive model for HTS.

# Fields
- `ec::Float64`: Critical electric field [V/m] (default: 1e-4)
- `jc::Float64`: Critical current density [A/m²]
- `n::Int`: Power-law exponent (higher = sharper transition)

# Theory
The power-law relates electric field to current density:

    **E** = E_c (|J|/J_c)^n (J/|J|)

which gives a nonlinear resistivity:

    ρ(|J|) = (E_c / J_c) * (|J| / J_c)^(n-1)

# References
- Rhyner (1993), Physica C, 212(3-4), 292-300.
"""
struct PowerLawMaterial <: AbstractHTSMaterial
    ec::Float64
    jc::Float64
    n::Int

    function PowerLawMaterial(; ec=E_C_DEFAULT, jc=JC_DEFAULT, n=N_DEFAULT)
        ec > 0 || throw(ArgumentError("E_c must be positive"))
        jc > 0 || throw(ArgumentError("J_c must be positive"))
        n > 0  || throw(ArgumentError("n must be positive"))
        new(ec, jc, n)
    end
end

"""
    PowerLawMaterial(jc, n; ec=E_C_DEFAULT)

Convenience constructor specifying J_c and n directly.
"""
PowerLawMaterial(jc::Real, n::Integer; ec=E_C_DEFAULT) =
    PowerLawMaterial(; ec=Float64(ec), jc=Float64(jc), n=Int(n))

"""
    resistivity(mat::PowerLawMaterial, J_norm::Real) -> Float64

Compute the nonlinear resistivity ρ(|J|) for the power-law model:

    ρ = (E_c / J_c) * (|J| / J_c)^(n-1)

# Arguments
- `mat::PowerLawMaterial`: Material model
- `J_norm::Real`: Magnitude of current density |J| [A/m²]

# Returns
- Resistivity ρ [Ω·m]

# Notes
Uses regularization (`EPS_REG`) to handle J_norm ≈ 0.
"""
function resistivity(mat::PowerLawMaterial, J_norm::Real)
    J_reg = max(J_norm, EPS_REG)
    return (mat.ec / mat.jc) * (J_reg / mat.jc)^(mat.n - 1)
end

"""
    electric_field(mat::PowerLawMaterial, J::VectorValue) -> VectorValue

Compute the electric field **E** from the power-law:

    **E** = E_c * (|J| / J_c)^n * (J / |J|)

# Arguments
- `mat::PowerLawMaterial`: Material model
- `J::VectorValue`: Current density vector [A/m²]

# Returns
- Electric field vector **E** [V/m]
"""
function electric_field(mat::PowerLawMaterial, J::VectorValue)
    J_norm = norm_safe(J)
    rho = resistivity(mat, J_norm)
    return rho * J
end

"""
    resistivity_derivative(mat::PowerLawMaterial, J_norm::Real) -> Float64

Compute dρ/d|J| for the power-law model (needed for Jacobian assembly):

    dρ/d|J| = (E_c / J_c²) * (n-1) * (|J| / J_c)^(n-2)

# Arguments
- `mat::PowerLawMaterial`: Material model
- `J_norm::Real`: Magnitude of current density |J| [A/m²]

# Returns
- Derivative of resistivity dρ/d|J| [Ω·m²/A]
"""
function resistivity_derivative(mat::PowerLawMaterial, J_norm::Real)
    J_reg = max(J_norm, EPS_REG)
    return (mat.ec / mat.jc^2) * (mat.n - 1) * (J_reg / mat.jc)^(mat.n - 2)
end
