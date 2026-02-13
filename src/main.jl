# ──────────────────────────────────────────────
# Main Driver Function
# ──────────────────────────────────────────────
#
# Following GridapMHD.jl pattern: a single `main(params)` entry point
# that dispatches to the appropriate formulation based on parameters.
# ──────────────────────────────────────────────

"""
    main(params::Dict; output::Dict=Dict{Symbol,Any}()) -> (xh, fullparams, output)

Main entry point for GridapHTS simulations.

Dispatches to the appropriate formulation solver based on
`params[:formulation]`:
- `:A`  -> A-formulation (magnetostatic vector potential)
- `:TA` -> T-A formulation (coupled HTS)

# Arguments
- `params::Dict`: Simulation parameters (see [`default_params`](@ref))
- `output::Dict`: Optional dictionary to collect solver statistics

# Returns
- `xh`: FE solution (scalar or MultiField)
- `fullparams::Dict`: Complete parameter dictionary with defaults filled in
- `output::Dict`: Solver statistics (iteration counts, residuals, timing)

# Example
```julia
params = Dict(
    :formulation => :A,
    :mesh => Dict(:type => :cartesian, :domain => (0,1,0,1), :partition => (20,20)),
    :fe_order => 1,
)
xh, fullparams, info = GridapHTS.main(params)
```
"""
function main(params::Dict; output::Dict=Dict{Symbol,Any}())
    # ── 1. Parameter processing ──────────────
    fullparams = add_default_params(params)
    validate_params(fullparams)

    formulation = fullparams[:formulation]
    @info "GridapHTS: Running $(formulation)-formulation"

    t_start = time()

    # ── 2. Dispatch to formulation ───────────
    xh, setup = _solve_formulation(Val(formulation), fullparams)

    t_elapsed = time() - t_start

    # ── 3. Post-processing ───────────────────
    output[:elapsed_time] = t_elapsed
    output[:formulation] = formulation

    # VTK output
    if get_nested(fullparams, :output, :vtk; default=false)
        _write_vtk(xh, setup, fullparams)
    end

    @info "GridapHTS: Completed in $(round(t_elapsed, digits=2))s"

    return xh, fullparams, output
end

# ──────────────────────────────────────────────
# Formulation dispatch via Val types
# ──────────────────────────────────────────────

"""Dispatch A-formulation."""
function _solve_formulation(::Val{:A}, params::Dict)
    uh, setup = solve_a_formulation(params)
    return uh, setup
end

"""Dispatch T-A formulation."""
function _solve_formulation(::Val{:TA}, params::Dict)
    xh, setup = solve_ta_formulation(params)
    return xh, setup
end

# ──────────────────────────────────────────────
# VTK Output
# ──────────────────────────────────────────────

"""Write solution to VTK file for visualization in ParaView."""
function _write_vtk(xh, setup, params)
    vtk_path = get_nested(params, :output, :vtk_path; default=nothing)

    if isnothing(vtk_path)
        problem_name = get(params, :problem_name, "gridaphts")
        vtk_path = problem_name
    end

    # Ensure output directory exists
    vtk_dir = dirname(vtk_path)
    if !isempty(vtk_dir)
        ensure_directory(vtk_dir)
    end

    Ω = setup.Ω
    formulation = params[:formulation]

    if formulation == :A
        writevtk(Ω, vtk_path, cellfields=["A" => xh])
        @info "VTK written to: $(vtk_path).vtu"
    elseif formulation == :TA
        A, T = xh
        writevtk(Ω, vtk_path, cellfields=["A" => A, "T" => T])
        @info "VTK written to: $(vtk_path).vtu"
    end
end
