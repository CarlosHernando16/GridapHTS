# ──────────────────────────────────────────────
# Newton-Raphson Solver with n-Continuation
# ──────────────────────────────────────────────
#
# Custom Newton-Raphson implementation tailored for HTS nonlinear problems.
# Includes n-continuation strategy: start with low n (smoother problem)
# and gradually increase to target n (sharper power-law).
# ──────────────────────────────────────────────

"""
    solve_newton(op; max_iter=50, rtol=1e-8, atol=1e-12, initial_guess=nothing)

Solve a nonlinear FE operator using Gridap's built-in Newton-Raphson.

# Arguments
- `op`: Nonlinear FE operator (from `FEOperator`)
- `max_iter::Int`: Maximum Newton iterations
- `rtol::Float64`: Relative tolerance on residual norm
- `atol::Float64`: Absolute tolerance on residual norm
- `initial_guess`: Optional FE function used to warm-start Newton

# Returns
- `uh`: FE solution

# Example
```julia
uh = solve_newton(setup.op; max_iter=50, rtol=1e-8)
```
"""
function solve_newton(op;
    max_iter::Int = 50,
    rtol::Float64 = 1e-8,
    atol::Float64 = 1e-12,
    initial_guess = nothing,
)
    nls = NLSolver(;
        show_trace = true,
        method     = :newton,
        iterations = max_iter,
        ftol       = atol,
        xtol       = rtol,
    )
    solver = FESolver(nls)

    if isnothing(initial_guess)
        uh = solve(solver, op)
    else
        # Warm-start nonlinear solve from the previous continuation step.
        uh = deepcopy(initial_guess)
        solve!(uh, solver, op)
    end

    return uh
end

"""
    solve_with_continuation(setup, params, n_schedule)

Solve a nonlinear HTS problem using n-continuation.

Strategy: start with a low power-law exponent (smooth problem), solve,
then use the solution as initial guess for the next higher n, repeating
until the target exponent is reached.

This significantly improves convergence for high-n problems (n > 20).

# Arguments
- `setup`: Named tuple from `setup_ta_formulation`
- `params::Dict`: Full parameter dictionary
- `n_schedule::Vector{Int}`: Sequence of n values, e.g., `[5, 10, 15, 20, 25]`

# Returns
- `xh`: Final FE solution at the target n

# Example
```julia
n_schedule = [5, 10, 15, 20, 25]
xh = solve_with_continuation(setup, params, n_schedule)
```

# Notes
- Each step rebuilds the nonlinear operator with updated material properties
- The previous solution is reused as the initial guess for the next step
"""
function solve_with_continuation(setup, params, n_schedule)
    @info "Starting n-continuation with schedule: $n_schedule"

    xh = nothing
    for (i, n) in enumerate(n_schedule)
        @info "  Step $i/$(length(n_schedule)): n = $n"

        # Update material parameters
        params_step = deepcopy(params)
        params_step[:material][:n_exponent] = n

        # Rebuild material and operator
        mat = material_from_params(params_step)
        μ_inv = setup.μ_inv
        D = num_cell_dims(setup.model)
        coordinate_system = get_nested(params_step, :mesh, :coordinate_system; default=:cartesian2d)

        res = ta_residual(mat, μ_inv, setup.dΩ, setup.dΩ_sc, D;
            coordinate_system=coordinate_system)
        jac = ta_jacobian(mat, μ_inv, setup.dΩ, setup.dΩ_sc, D;
            coordinate_system=coordinate_system)
        op = FEOperator(res, jac, setup.X, setup.Y)

        # Solve
        solver_params = params[:solver]
        xh = solve_newton(op;
            max_iter = solver_params[:max_iter],
            rtol     = solver_params[:rtol],
            atol     = get(solver_params, :atol, 1e-12),
            initial_guess = xh,
        )

        @info "  Step $i complete."
    end

    return xh
end
