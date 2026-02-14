#!/usr/bin/env julia

# ──────────────────────────────────────────────
# Phase 1 Validation: 2D Axisymmetric Coil (A-Formulation)
# ──────────────────────────────────────────────
#
# Solves the magnetostatic A-formulation on an (r,z) domain using the
# axisymmetric 2D mode:
#   (1/r) ∂_r ( r (1/μ) ∂_r A ) + ∂_z( (1/μ) ∂_z A ) = -J_phi
#
# GridapHTS handles the axisymmetric weighting (2πr) internally when
# :mesh[:coordinate_system] == :axisymmetric2d.
#
# This example uses a simple rectangular annular coil region with uniform
# azimuthal source current density J_phi derived from a target current I_coil.
# ──────────────────────────────────────────────

using DrWatson
@quickactivate "GridapHTS"

using GridapHTS

# ── Coil & domain parameters ──────────────────
I_coil = 1.0e3  # [A] total current assigned to the coil cross-section

# Axisymmetric (r,z) geometry [m]
r_inner = 0.020
r_outer = 0.030
z_min = -0.005
z_max = 0.005

# Computational air box around coil [m]
r_domain_max = 0.080
z_domain_half = 0.050

# Uniform J_phi in the coil cross-section
coil_area = (r_outer - r_inner) * (z_max - z_min)
J_coil = I_coil / coil_area  # [A/m^2]

coil_source(x) = begin
    r = x[1]
    z = x[2]
    in_coil = (r_inner <= r <= r_outer) && (z_min <= z <= z_max)
    return in_coil ? J_coil : 0.0
end

# ── Solver parameters ─────────────────────────
params = Dict(
    :problem_name => "axisymmetric_coil_2d",
    :formulation  => :A,
    :fe_order     => 1,

    :mesh => Dict(
        :type              => :cartesian,
        :coordinate_system => :axisymmetric2d,
        :domain            => (0.0, r_domain_max, -z_domain_half, z_domain_half),
        :partition         => (120, 120),
    ),

    :material => Dict(
        :mu => GridapHTS.MU_0,
    ),

    :bcs => Dict(
        :dirichlet_tags => ["boundary"],
        :dirichlet_values => 0.0,
    ),

    :source => coil_source,

    :output => Dict(
        :vtk      => true,
        :vtk_path => joinpath(plotsdir("phase1"), "axisymmetric_coil_2d"),
    ),
)

# ── Run simulation ────────────────────────────
@info "Running axisymmetric 2D coil example..."
xh, fullparams, info = GridapHTS.main(params)

# ── Save metadata / derived quantities ────────
results = Dict(
    "problem_name"      => fullparams[:problem_name],
    "formulation"       => String(fullparams[:formulation]),
    "coordinate_system" => String(fullparams[:mesh][:coordinate_system]),
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
    "mesh_partition"    => fullparams[:mesh][:partition],
)

result_path = datadir("phase1", "axisymmetric_coil_2d_results.jld2")
mkpath(dirname(result_path))
@tagsave(result_path, results)

@info "Axisymmetric 2D coil run finished."
@info "  Elapsed time: $(round(info[:elapsed_time], digits=3)) s"
@info "  VTK: $(fullparams[:output][:vtk_path]).vtu"
@info "  JLD2: $result_path"
