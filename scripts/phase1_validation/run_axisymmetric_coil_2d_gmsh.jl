#!/usr/bin/env julia

# ──────────────────────────────────────────────
# Phase 1 Validation: 2D Axisymmetric Coil (Gmsh-tagged source)
# ──────────────────────────────────────────────
#
# This variant assigns source current by mesh subdomain tag instead of
# coordinate windows. The mesh must include:
#   - Physical surface tag: "coil"      (current-carrying region)
#   - Physical boundary tag: "boundary" (Dirichlet boundary for A)
#
# The solver integrates J only over the "coil" subdomain through
# `:source_tag => "coil"`.
# ──────────────────────────────────────────────

using DrWatson
@quickactivate "GridapHTS"

using GridapHTS

# ── Inputs ────────────────────────────────────
I_coil = 1.0e3 # [A]

# Geometric metadata used only to compute a uniform J in the tagged coil.
# Keep these consistent with your Gmsh geometry for physical correctness.
r_inner = 0.020
r_outer = 0.030
z_min = -0.005
z_max = 0.005
coil_area = (r_outer - r_inner) * (z_max - z_min)
J_coil = I_coil / coil_area # [A/m^2]

# Default mesh path (override with ENV if desired)
mesh_file = get(
    ENV,
    "GRIDAPHTS_AXI_GMSH_FILE",
    datadir("meshes", "benchmarks", "axisymmetric_coil_2d_tagged.msh"),
)

isfile(mesh_file) || error(
    "Mesh file not found: $mesh_file\n" *
    "Provide a Gmsh mesh with Physical Surface(\"coil\") and " *
    "Physical Curve(\"boundary\").\n" *
    "You can override the path with ENV[\"GRIDAPHTS_AXI_GMSH_FILE\"]."
)

# ── Solver setup ──────────────────────────────
params = Dict(
    :problem_name => "axisymmetric_coil_2d_gmsh",
    :formulation  => :A,
    :fe_order     => 1,

    :mesh => Dict(
        :type              => :gmsh,
        :file              => mesh_file,
        :coordinate_system => :axisymmetric2d,
    ),

    :material => Dict(
        :mu => GridapHTS.MU_0,
    ),

    :bcs => Dict(
        :dirichlet_tags   => ["boundary"],
        :dirichlet_values => 0.0,
    ),

    # Source applied only on `source_tag` subdomain.
    :source     => J_coil,
    :source_tag => "coil",

    :output => Dict(
        :vtk      => true,
        :vtk_path => joinpath(plotsdir("phase1"), "axisymmetric_coil_2d_gmsh"),
    ),
)

# ── Run ───────────────────────────────────────
@info "Running axisymmetric 2D coil (Gmsh-tagged source) example..."
xh, fullparams, info = GridapHTS.main(params)

# ── Save metadata ─────────────────────────────
results = Dict(
    "problem_name"      => fullparams[:problem_name],
    "formulation"       => String(fullparams[:formulation]),
    "coordinate_system" => String(fullparams[:mesh][:coordinate_system]),
    "mesh_file"         => mesh_file,
    "source_tag"        => fullparams[:source_tag],
    "elapsed_time_s"    => info[:elapsed_time],
    "I_coil_A"          => I_coil,
    "J_coil_A_per_m2"   => J_coil,
    "coil_area_m2"      => coil_area,
    "coil_region"       => Dict(
        "r_inner_m" => r_inner,
        "r_outer_m" => r_outer,
        "z_min_m"   => z_min,
        "z_max_m"   => z_max,
    ),
)

result_path = datadir("phase1", "axisymmetric_coil_2d_gmsh_results.jld2")
mkpath(dirname(result_path))
@tagsave(result_path, results)

@info "Axisymmetric 2D Gmsh-tagged run finished."
@info "  Elapsed time: $(round(info[:elapsed_time], digits=3)) s"
@info "  VTK: $(fullparams[:output][:vtk_path]).vtu"
@info "  JLD2: $result_path"
