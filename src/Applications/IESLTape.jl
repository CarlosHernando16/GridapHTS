# ──────────────────────────────────────────────
# IESL 2D Rectangular Tape Benchmark
# ──────────────────────────────────────────────
#
# Standard benchmark from htsmodelling.com:
# 2D rectangular HTS tape in an applied magnetic field.
# Validates screening current distribution against reference data.
# ──────────────────────────────────────────────

"""
    setup_iesl_benchmark(; kwargs...) -> Dict{Symbol,Any}

Create parameter dictionary for the IESL 2D rectangular tape benchmark.

This benchmark models a single HTS tape cross-section under an applied
magnetic field, computing the screening current distribution and
comparing against reference data from htsmodelling.com.

# Keyword Arguments
- `mesh_file::Union{String,Nothing}`: Path to .msh file (default: uses Cartesian)
- `b_max::Float64`: Maximum applied field [T] (default: 1.0)
- `jc::Float64`: Critical current density [A/m²] (default: 3e10)
- `n::Int`: Power-law exponent (default: 21)
- `fe_order::Int`: FE polynomial order (default: 1)
- `partition::Tuple`: Mesh partition for Cartesian mesh (default: (40, 10))

# Returns
- `params::Dict{Symbol,Any}`: Complete parameter dictionary ready for `main()`

# Example
```julia
params = setup_iesl_benchmark(; jc=3e10, n=21, b_max=0.5)
xh, fullparams, info = GridapHTS.main(params)
```

# References
- IESL HTS Modelling Workgroup: https://htsmodelling.com/
- Benchmark Problem 1: 2D Tape in Applied Field
"""
function setup_iesl_benchmark(;
    mesh_file::Union{String,Nothing} = nothing,
    b_max::Float64 = 1.0,
    jc::Float64 = 3e10,
    n::Int = 21,
    fe_order::Int = 1,
    partition::Tuple = (40, 10),
)
    # Tape geometry: width = 12 mm, thickness = 1 μm
    # Modeled domain includes surrounding air region
    tape_width = 12e-3      # [m]
    tape_thickness = 1e-6   # [m]

    # Domain: air region around tape (10x tape dimensions)
    air_margin = 10 * tape_width
    domain = (
        -air_margin, air_margin,
        -air_margin, air_margin,
    )

    mesh_params = if isnothing(mesh_file)
        Dict{Symbol,Any}(
            :type      => :cartesian,
            :domain    => domain,
            :partition => partition,
        )
    else
        Dict{Symbol,Any}(
            :type => :gmsh,
            :file => mesh_file,
            :superconductor_tag => "tape",
        )
    end

    params = Dict{Symbol,Any}(
        :problem_name => "IESL_2D_Tape",
        :formulation  => :TA,
        :fe_order     => fe_order,

        :mesh => mesh_params,

        :material => Dict{Symbol,Any}(
            :type       => :power_law,
            :jc         => jc,
            :n_exponent => n,
            :ec         => E_C_DEFAULT,
            :mu         => MU_0,
        ),

        :bcs => Dict{Symbol,Any}(
            :dirichlet_tags   => ["boundary"],
            :dirichlet_values => nothing,
            :b_applied        => b_max,
        ),

        :solver => Dict{Symbol,Any}(
            :type           => :newton,
            :max_iter       => 50,
            :rtol           => 1e-6,
            :atol           => 1e-12,
            :n_continuation => [5, 10, 15, n],
        ),

        :output => Dict{Symbol,Any}(
            :vtk            => true,
            :vtk_path       => nothing,
            :save_residuals => true,
        ),

        :gauge => Dict{Symbol,Any}(
            :type => :coulomb,
        ),

        # IESL-specific metadata
        :benchmark => Dict{Symbol,Any}(
            :name           => "IESL 2D Rectangular Tape",
            :tape_width     => tape_width,
            :tape_thickness => tape_thickness,
            :reference_url  => "https://htsmodelling.com/",
        ),
    )

    return params
end
