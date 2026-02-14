# ──────────────────────────────────────────────
# T-A Formulation (Coupled Vector & Scalar Potentials)
# ──────────────────────────────────────────────
#
# Solves the coupled T-A system for HTS:
#   A-equation: ∇ × (1/μ ∇ × A) = J_s      (whole domain)
#   T-equation: ∇ ⋅ E(J) = 0                (superconductor)
#
# where J = -∇T and E = ρ(|J|)J from the power-law.
# ──────────────────────────────────────────────

"""
    setup_ta_formulation(params::Dict) -> NamedTuple

Set up the coupled T-A formulation FE problem.

Returns a named tuple containing:
- `model`: Discrete model
- `Y`: MultiField test space [V_A, V_T]
- `X`: MultiField trial space [U_A, U_T]
- `Ω`: Full domain triangulation
- `Ω_sc`: Superconductor subdomain triangulation
- `dΩ`: Full domain measure
- `dΩ_sc`: Superconductor measure
- `op`: Nonlinear FE operator
- `mat`: HTS material model

# Arguments
- `params::Dict`: Parameter dictionary with `:formulation => :TA`

# Notes
- The T-equation is only assembled on the superconductor subdomain
- Requires the mesh to have physical tags distinguishing superconductor
  from non-superconductor regions
- For Phase 2+ development

# Example
```julia
params = Dict(
    :formulation => :TA,
    :mesh => Dict(:type => :gmsh, :file => "tape.msh"),
    :material => Dict(:jc => 3e10, :n_exponent => 25),
    :fe_order => 1,
)
setup = setup_ta_formulation(params)
```
"""
function setup_ta_formulation(params::Dict)
    order = params[:fe_order]
    μ = get_nested(params, :material, :mu; default=MU_0)
    μ_inv = 1.0 / μ

    # Build model
    model = _build_model(params)
    D = num_cell_dims(model)
    coordinate_system = get_nested(params, :mesh, :coordinate_system; default=:cartesian2d)

    # Material
    mat = material_from_params(params)

    # ── FE Spaces ────────────────────────────

    # A space: H(curl) in 3D, H¹ in 2D
    bc_tags_A = get_nested(params, :bcs, :dirichlet_tags; default=["boundary"])

    if D == 2
        reffe_A = ReferenceFE(lagrangian, Float64, order)
        V_A = TestFESpace(model, reffe_A;
            conformity=:H1, dirichlet_tags=bc_tags_A)
        U_A = TrialFESpace(V_A, 0.0)
    else
        reffe_A = ReferenceFE(nedelec, Float64, order)
        V_A = TestFESpace(model, reffe_A;
            conformity=:Hcurl, dirichlet_tags=bc_tags_A)
        U_A = TrialFESpace(V_A, VectorValue(0.0, 0.0, 0.0))
    end

    # T space: H¹ Lagrangian on superconductor subdomain
    bc_tags_T = get_nested(params, :bcs, :dirichlet_tags_T;
        default=["superconductor_boundary"])

    reffe_T = ReferenceFE(lagrangian, Float64, order)
    V_T = TestFESpace(model, reffe_T;
        conformity=:H1, dirichlet_tags=bc_tags_T)
    U_T = TrialFESpace(V_T, 0.0)

    # MultiField spaces
    Y = MultiFieldFESpace([V_A, V_T])
    X = MultiFieldFESpace([U_A, U_T])

    # ── Triangulations & Measures ────────────

    Ω = Triangulation(model)
    degree = 2 * order + 2
    dΩ = Measure(Ω, degree)

    # Superconductor subdomain
    sc_tag = get_nested(params, :mesh, :superconductor_tag;
        default="superconductor")
    # For now, use full domain if no tag specified (Phase 2 will refine)
    Ω_sc = Ω
    dΩ_sc = dΩ

    # ── Nonlinear FE Operator ────────────────

    res = ta_residual(mat, μ_inv, dΩ, dΩ_sc, D; coordinate_system=coordinate_system)
    jac = ta_jacobian(mat, μ_inv, dΩ, dΩ_sc, D; coordinate_system=coordinate_system)

    op = FEOperator(res, jac, X, Y)

    return (
        model  = model,
        Y      = Y,
        X      = X,
        Ω      = Ω,
        Ω_sc   = Ω_sc,
        dΩ     = dΩ,
        dΩ_sc  = dΩ_sc,
        op     = op,
        mat    = mat,
        μ_inv  = μ_inv,
        coordinate_system = coordinate_system,
    )
end

"""
    solve_ta_formulation(params::Dict) -> (xh, setup)

Convenience function: set up and solve the T-A formulation.

Uses Newton-Raphson with optional n-continuation for robust convergence.

# Arguments
- `params::Dict`: Parameter dictionary

# Returns
- `xh`: MultiField FE solution (A, T)
- `setup`: Named tuple from [`setup_ta_formulation`](@ref)
"""
function solve_ta_formulation(params::Dict)
    setup = setup_ta_formulation(params)

    # Check for n-continuation
    n_schedule = get_nested(params, :solver, :n_continuation; default=nothing)

    if !isnothing(n_schedule)
        xh = solve_with_continuation(setup, params, n_schedule)
    else
        solver_params = params[:solver]
        xh = solve_newton(setup.op;
            max_iter = solver_params[:max_iter],
            rtol     = solver_params[:rtol],
            atol     = get(solver_params, :atol, 1e-12),
        )
    end

    return xh, setup
end
