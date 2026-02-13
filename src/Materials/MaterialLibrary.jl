# ──────────────────────────────────────────────
# Material Library: Preset HTS Material Definitions
# ──────────────────────────────────────────────

"""
    rebco_default(; kwargs...) -> PowerLawMaterial

Default REBCO (REBa₂Cu₃O₇₋ₓ) tape material parameters.

Typical values at 77 K, self-field:
- J_c = 3×10¹⁰ A/m² (corresponds to ~300 A/cm-width for 1 μm film)
- n = 25
- E_c = 10⁻⁴ V/m

# Keyword Arguments
Override any default:
- `jc::Real`: Critical current density [A/m²]
- `n::Integer`: Power-law exponent
- `ec::Real`: Critical electric field [V/m]
"""
function rebco_default(; jc=3e10, n=25, ec=E_C_DEFAULT)
    return PowerLawMaterial(; ec=Float64(ec), jc=Float64(jc), n=Int(n))
end

"""
    bscco_default(; kwargs...) -> PowerLawMaterial

Default Bi-2223 (BiSCCO) material parameters.

Typical values at 77 K, self-field:
- J_c = 1×10⁹ A/m²
- n = 15
- E_c = 10⁻⁴ V/m

# Keyword Arguments
Override any default:
- `jc::Real`: Critical current density [A/m²]
- `n::Integer`: Power-law exponent
- `ec::Real`: Critical electric field [V/m]
"""
function bscco_default(; jc=1e9, n=15, ec=E_C_DEFAULT)
    return PowerLawMaterial(; ec=Float64(ec), jc=Float64(jc), n=Int(n))
end

"""
    material_from_params(params::Dict) -> AbstractHTSMaterial

Construct an HTS material from a parameter dictionary.

Dispatches on `params[:material][:type]`:
- `:power_law` -> `PowerLawMaterial`
- `:field_dependent` -> `FieldDependentMaterial`

# Example
```julia
params = Dict(:material => Dict(:type => :power_law, :jc => 1e9, :n_exponent => 25))
mat = material_from_params(params)
```
"""
function material_from_params(params::Dict)
    mp = params[:material]
    mat_type = get(mp, :type, :power_law)

    if mat_type == :power_law
        return PowerLawMaterial(;
            ec = get(mp, :ec, E_C_DEFAULT),
            jc = get(mp, :jc, JC_DEFAULT),
            n  = get(mp, :n_exponent, N_DEFAULT),
        )
    elseif mat_type == :field_dependent
        jc_model = KimModel(;
            jc0 = get(mp, :jc, JC_DEFAULT),
            b0  = get(mp, :b0, 0.1),
        )
        return FieldDependentMaterial(;
            ec       = get(mp, :ec, E_C_DEFAULT),
            n        = get(mp, :n_exponent, N_DEFAULT),
            jc_model = jc_model,
        )
    else
        throw(ArgumentError("Unknown material type: $mat_type"))
    end
end
