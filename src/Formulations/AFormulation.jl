# ──────────────────────────────────────────────
# A-Formulation (Magnetic Vector Potential)
# ──────────────────────────────────────────────
#
# Solves the magnetostatic problem:
#   ∇ × (1/μ ∇ × A) = J_s
#
# with appropriate boundary conditions and gauge condition.
# ──────────────────────────────────────────────

"""
    setup_a_formulation(params::Dict) -> NamedTuple

Set up the A-formulation FE problem from a parameter dictionary.

Returns a named tuple containing:
- `model`: Discrete model
- `V`: Test FE space (H(curl) or H¹ depending on dimension)
- `U`: Trial FE space
- `Ω`: Triangulation
- `dΩ`: Measure
- `op`: FE operator

# Arguments
- `params::Dict`: Parameter dictionary (see [`default_params`](@ref))

# Example
```julia
params = Dict(
    :formulation => :A,
    :mesh => Dict(:type => :cartesian, :domain => (0,1,0,1), :partition => (20,20)),
    :material => Dict(:mu => 4π*1e-7),
    :bcs => Dict(:dirichlet_tags => ["boundary"]),
    :fe_order => 1,
)
setup = setup_a_formulation(params)
```
"""
function setup_a_formulation(params::Dict)
    # Extract parameters
    order = params[:fe_order]
    μ = get_nested(params, :material, :mu; default=MU_0)
    μ_inv = 1.0 / μ

    # Build discrete model
    model = _build_model(params)
    coordinate_system = get_nested(params, :mesh, :coordinate_system; default=:cartesian2d)

    # FE spaces
    # For 2D scalar A (out-of-plane component): use Lagrangian H¹
    # For 3D vector A: use Nédélec H(curl)
    D = num_cell_dims(model)

    if D == 2
        reffe = ReferenceFE(lagrangian, Float64, order)
        conformity = :H1
    else
        reffe = ReferenceFE(nedelec, Float64, order)
        conformity = :Hcurl
    end

    bc_tags = get_nested(params, :bcs, :dirichlet_tags; default=["boundary"])
    bc_vals = get_nested(params, :bcs, :dirichlet_values; default=nothing)

    if isnothing(bc_vals)
        bc_vals = D == 2 ? 0.0 : VectorValue(0.0, 0.0, 0.0)
    end

    V = TestFESpace(model, reffe; conformity=conformity, dirichlet_tags=bc_tags)
    U = TrialFESpace(V, bc_vals)

    # Triangulation and measures
    Ω = Triangulation(model)
    degree = 2 * order + 1
    dΩ = Measure(Ω, degree)

    # Source term (zero by default)
    source = get_nested(params, :source; default=nothing)
    if isnothing(source)
        source = D == 2 ? (x -> 0.0) : (x -> VectorValue(0.0, 0.0, 0.0))
    end
    source_tag = get_nested(params, :source_tag; default=nothing)

    # Weak forms
    a = a_bilinear_form(μ_inv, dΩ, D; coordinate_system=coordinate_system)
    if isnothing(source_tag)
        l = a_linear_form(source, dΩ, D; coordinate_system=coordinate_system)
    else
        Ω_source = Triangulation(model, source_tag)
        dΩ_source = Measure(Ω_source, degree)
        source_fun = source isa Number ? (x -> source) : source
        l = a_linear_form(source_fun, dΩ_source, D; coordinate_system=coordinate_system)
    end

    # FE Operator
    op = AffineFEOperator(a, l, U, V)

    return (
        model = model,
        V = V,
        U = U,
        Ω = Ω,
        dΩ = dΩ,
        op = op,
        μ_inv = μ_inv,
        coordinate_system = coordinate_system,
        source_tag = source_tag,
    )
end

"""
    solve_a_formulation(params::Dict) -> (uh, setup)

Convenience function: set up and solve the A-formulation in one call.

# Arguments
- `params::Dict`: Parameter dictionary

# Returns
- `uh`: FE solution (magnetic vector potential A)
- `setup`: Named tuple from [`setup_a_formulation`](@ref)
"""
function solve_a_formulation(params::Dict)
    setup = setup_a_formulation(params)
    uh = solve(setup.op)
    return uh, setup
end

# ── Internal helpers ─────────────────────────

"""Build a DiscreteModel from the mesh parameters."""
function _build_model(params::Dict)
    mesh = params[:mesh]
    if mesh[:type] == :cartesian
        domain = mesh[:domain]
        partition = mesh[:partition]
        return CartesianDiscreteModel(domain, partition)
    elseif mesh[:type] == :gmsh
        return load_gmsh_model(mesh[:file])
    else
        error("Unknown mesh type: $(mesh[:type])")
    end
end
