# ──────────────────────────────────────────────
# Phase 1 Validation: A-Formulation (Magnetostatic)
# ──────────────────────────────────────────────
#
# Solves the 2D magnetostatic problem using the A-formulation:
#   ∇ × (1/μ ∇ × A) = J_s
#
# In 2D, A is a scalar (out-of-plane component) and this reduces to:
#   -∇ · (1/μ ∇A) = J_s
#
# Test case: prescribed source current, compare B = ∇ × A field.
# ──────────────────────────────────────────────

using DrWatson
@quickactivate "GridapHTS"

using GridapHTS
using Gridap

# ── Parameters ───────────────────────────────
params = Dict(
    :formulation => :A,
    :fe_order    => 1,

    :mesh => Dict(
        :type      => :cartesian,
        :domain    => (0.0, 1.0, 0.0, 1.0),
        :partition => (40, 40),
    ),

    :material => Dict(
        :mu => GridapHTS.MU_0,
    ),

    :bcs => Dict(
        :dirichlet_tags => ["boundary"],
    ),

    :output => Dict(
        :vtk      => true,
        :vtk_path => joinpath(plotsdir("phase1"), "a_formulation"),
    ),
)

# ── Solve ────────────────────────────────────
@info "Running A-formulation validation..."
xh, fullparams, info = GridapHTS.main(params)

@info "A-formulation completed"
@info "  Elapsed time: $(round(info[:elapsed_time], digits=3))s"

# ── Post-process ─────────────────────────────
# The solution is the magnetic vector potential A (scalar in 2D)
# B = curl(A) = (∂A/∂y, -∂A/∂x) in 2D

# Save results
results = Dict(
    "elapsed_time" => info[:elapsed_time],
    "formulation"  => "A",
    "fe_order"     => fullparams[:fe_order],
    "partition"    => fullparams[:mesh][:partition],
)
@tagsave(datadir("phase1", "a_formulation_results.jld2"), results)
@info "Results saved."
