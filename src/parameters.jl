# ──────────────────────────────────────────────
# Parameter Management
# ──────────────────────────────────────────────
#
# Following GridapMHD.jl pattern: all simulation configuration lives in a
# single Dict{Symbol,Any}. `add_default_params` fills in missing entries
# with sensible defaults so that scripts only need to specify overrides.
# ──────────────────────────────────────────────

"""
    default_params() -> Dict{Symbol,Any}

Return the complete default parameter dictionary for a GridapHTS simulation.
"""
function default_params()
    Dict{Symbol,Any}(
        # ── Formulation ──────────────────────────
        :formulation => :A,        # :A or :TA
        :fe_order    => 1,         # Polynomial order

        # ── Mesh ─────────────────────────────────
        :mesh => Dict{Symbol,Any}(
            :type      => :cartesian,  # :cartesian or :gmsh
            :coordinate_system => :cartesian2d, # :cartesian2d or :axisymmetric2d (for 2D)
            :domain    => (0, 1, 0, 1),
            :partition => (20, 20),
            :file      => nothing,     # Path to .msh file (for :gmsh)
        ),

        # ── Material ─────────────────────────────
        :material => Dict{Symbol,Any}(
            :type       => :power_law,
            :jc         => JC_DEFAULT,
            :n_exponent => N_DEFAULT,
            :ec         => E_C_DEFAULT,
            :mu         => MU_0,
        ),

        # ── Boundary Conditions ──────────────────
        :bcs => Dict{Symbol,Any}(
            :dirichlet_tags   => ["boundary"],
            :dirichlet_values => nothing,  # nothing => zero BC
            :neumann_tags     => String[],
            :neumann_values   => nothing,
        ),

        # ── Solver ───────────────────────────────
        :solver => Dict{Symbol,Any}(
            :type     => :newton,     # :newton, :picard
            :max_iter => 50,
            :rtol     => 1e-8,
            :atol     => 1e-12,
            :n_continuation => nothing,  # e.g., [5, 10, 15, 20, 25]
        ),

        # ── Output ───────────────────────────────
        :output => Dict{Symbol,Any}(
            :vtk        => false,
            :vtk_path   => nothing,  # Auto-generated if nothing
            :save_residuals => false,
        ),

        # ── Gauge ────────────────────────────────
        :gauge => Dict{Symbol,Any}(
            :type => :none,  # :none, :coulomb, :tree_cotree
        ),
    )
end

"""
    add_default_params(params::Dict) -> Dict{Symbol,Any}

Merge user-provided `params` with [`default_params()`](@ref), filling in
any missing entries with defaults.

# Arguments
- `params::Dict`: User-supplied parameter dictionary (may be partial)

# Returns
- Complete parameter dictionary with all required fields

# Example
```julia
user_params = Dict(:formulation => :TA, :material => Dict(:jc => 3e10))
full_params = add_default_params(user_params)
```
"""
function add_default_params(params::Dict)
    # Convert String keys to Symbol keys if needed
    sp = _symbolize_keys(params)
    return merge_nested(default_params(), sp)
end

"""
    validate_params(params::Dict)

Check that the parameter dictionary has internally consistent values.
Throws `ArgumentError` on invalid combinations.
"""
function validate_params(params::Dict)
    form = params[:formulation]
    form in (:A, :TA) || throw(ArgumentError(
        "Unknown formulation: $form. Expected :A or :TA"))

    mesh = params[:mesh]
    mesh[:type] in (:cartesian, :gmsh) || throw(ArgumentError(
        "Unknown mesh type: $(mesh[:type]). Expected :cartesian or :gmsh"))
    mesh[:coordinate_system] in (:cartesian2d, :axisymmetric2d) || throw(ArgumentError(
        "Unknown coordinate_system: $(mesh[:coordinate_system]). " *
        "Expected :cartesian2d or :axisymmetric2d"))

    if mesh[:type] == :gmsh
        isnothing(mesh[:file]) && throw(ArgumentError(
            "Mesh type :gmsh requires :file parameter"))
    end

    # Axisymmetric 2D requires radial coordinate r >= 0
    if mesh[:type] == :cartesian && mesh[:coordinate_system] == :axisymmetric2d
        domain = mesh[:domain]
        r_min = domain[1]
        r_min >= 0 || throw(ArgumentError(
            "Axisymmetric 2D requires radial coordinate r >= 0, got r_min=$r_min"))
    end

    mat = params[:material]
    mat[:jc] > 0 || throw(ArgumentError("J_c must be positive, got $(mat[:jc])"))
    mat[:n_exponent] > 0 || throw(ArgumentError(
        "n_exponent must be positive, got $(mat[:n_exponent])"))
    mat[:ec] > 0 || throw(ArgumentError("E_c must be positive, got $(mat[:ec])"))
    mat[:mu] > 0 || throw(ArgumentError("μ must be positive, got $(mat[:mu])"))

    solver = params[:solver]
    solver[:type] in (:newton, :picard) || throw(ArgumentError(
        "Unknown solver type: $(solver[:type])"))
    solver[:max_iter] > 0 || throw(ArgumentError("max_iter must be positive"))
    solver[:rtol] > 0 || throw(ArgumentError("rtol must be positive"))

    return nothing
end

# ── Internal helpers ─────────────────────────

"""Convert all Dict keys to Symbols recursively."""
function _symbolize_keys(d::Dict)
    result = Dict{Symbol,Any}()
    for (k, v) in d
        sk = k isa Symbol ? k : Symbol(k)
        result[sk] = v isa Dict ? _symbolize_keys(v) : v
    end
    return result
end
