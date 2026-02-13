# ──────────────────────────────────────────────
# Field-Dependent Critical Current Models
# ──────────────────────────────────────────────
#
# Models for J_c(B, T) dependence, essential for accurate HTS modeling
# under external magnetic fields and varying temperatures.
# ──────────────────────────────────────────────

"""
    KimModel

Kim model for field-dependent critical current density:

    J_c(B) = J_c0 / (1 + |B| / B_0)

# Fields
- `jc0::Float64`: Self-field critical current density [A/m²]
- `b0::Float64`: Characteristic field [T]

# References
- Kim et al. (1962), Phys. Rev. Lett., 9(7), 306.
"""
struct KimModel
    jc0::Float64
    b0::Float64

    function KimModel(; jc0=JC_DEFAULT, b0=0.1)
        jc0 > 0 || throw(ArgumentError("J_c0 must be positive"))
        b0 > 0  || throw(ArgumentError("B_0 must be positive"))
        new(jc0, b0)
    end
end

"""
    critical_current_density(model::KimModel, B::VectorValue) -> Float64

Evaluate J_c(B) using the Kim model.

# Arguments
- `model::KimModel`: Kim model parameters
- `B::VectorValue`: Magnetic flux density vector [T]

# Returns
- Critical current density J_c [A/m²]
"""
function critical_current_density(model::KimModel, B::VectorValue)
    B_norm = norm_safe(B)
    return model.jc0 / (1.0 + B_norm / model.b0)
end

"""
    critical_current_density(model::KimModel, B_norm::Real) -> Float64

Scalar version accepting |B| directly.
"""
function critical_current_density(model::KimModel, B_norm::Real)
    return model.jc0 / (1.0 + abs(B_norm) / model.b0)
end

"""
    FieldDependentMaterial <: AbstractHTSMaterial

Power-law material with field-dependent J_c(B).

Combines the standard power-law E-J relation with a field-dependent
critical current model:

    E = E_c * (|J| / J_c(B))^n * (J / |J|)

# Fields
- `ec::Float64`: Critical electric field [V/m]
- `n::Int`: Power-law exponent
- `jc_model::KimModel`: Field-dependent J_c model
"""
struct FieldDependentMaterial <: AbstractHTSMaterial
    ec::Float64
    n::Int
    jc_model::KimModel

    function FieldDependentMaterial(;
        ec=E_C_DEFAULT,
        n=N_DEFAULT,
        jc_model=KimModel()
    )
        ec > 0 || throw(ArgumentError("E_c must be positive"))
        n > 0  || throw(ArgumentError("n must be positive"))
        new(ec, n, jc_model)
    end
end

"""
    resistivity(mat::FieldDependentMaterial, J_norm::Real, B::VectorValue) -> Float64

Compute the nonlinear resistivity with field-dependent J_c:

    ρ = (E_c / J_c(B)) * (|J| / J_c(B))^(n-1)

# Arguments
- `mat::FieldDependentMaterial`: Material model
- `J_norm::Real`: Magnitude of current density |J| [A/m²]
- `B::VectorValue`: Magnetic flux density vector [T]
"""
function resistivity(mat::FieldDependentMaterial, J_norm::Real, B::VectorValue)
    jc = critical_current_density(mat.jc_model, B)
    J_reg = max(J_norm, EPS_REG)
    return (mat.ec / jc) * (J_reg / jc)^(mat.n - 1)
end

"""
    electric_field(mat::FieldDependentMaterial, J::VectorValue, B::VectorValue) -> VectorValue

Compute the electric field with field-dependent J_c.
"""
function electric_field(mat::FieldDependentMaterial, J::VectorValue, B::VectorValue)
    J_norm = norm_safe(J)
    rho = resistivity(mat, J_norm, B)
    return rho * J
end
