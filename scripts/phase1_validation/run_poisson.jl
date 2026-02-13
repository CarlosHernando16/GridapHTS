# ──────────────────────────────────────────────
# Phase 1 Validation: Poisson Equation
# ──────────────────────────────────────────────
#
# Validates the Gridap setup by solving a simple Poisson problem:
#   -Δu = f on [0,1]²
#   u = 0 on ∂Ω
#
# Analytical solution: u(x,y) = sin(πx) sin(πy)
# Source: f(x,y) = 2π² sin(πx) sin(πy)
# ──────────────────────────────────────────────

using DrWatson
@quickactivate "GridapHTS"

using Gridap

# ── Parameters ───────────────────────────────
domain = (0, 1, 0, 1)
partition = (20, 20)
order = 1

# ── Model & FE Space ────────────────────────
model = CartesianDiscreteModel(domain, partition)
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe; conformity=:H1, dirichlet_tags="boundary")
U = TrialFESpace(V, 0.0)

# ── Triangulation & Measures ─────────────────
Ω = Triangulation(model)
dΩ = Measure(Ω, 2 * order + 1)

# ── Weak Form ───────────────────────────────
a(u, v) = ∫(∇(u) ⋅ ∇(v))dΩ
f(x) = 2π^2 * sin(π * x[1]) * sin(π * x[2])
l(v) = ∫(f * v)dΩ

# ── Solve ────────────────────────────────────
op = AffineFEOperator(a, l, U, V)
uh = solve(op)

# ── Error Analysis ───────────────────────────
u_exact(x) = sin(π * x[1]) * sin(π * x[2])
e = u_exact - uh
l2_error = sqrt(sum(∫(e * e)dΩ))
h1_error = sqrt(sum(∫(∇(e) ⋅ ∇(e))dΩ))

@info "Poisson Validation Results"
@info "  Mesh: $(partition[1])×$(partition[2]), order=$order"
@info "  L² error: $l2_error"
@info "  H¹ error: $h1_error"

# ── VTK Output ───────────────────────────────
output_dir = plotsdir("phase1")
mkpath(output_dir)
output_file = joinpath(output_dir, "poisson")
writevtk(Ω, output_file, cellfields=["uh" => uh, "error" => e])
@info "VTK saved to: $(output_file).vtu"

# ── Save results with DrWatson ───────────────
results = @strdict l2_error h1_error partition order
@tagsave(datadir("phase1", "poisson_results.jld2"), results)
@info "Results saved to: $(datadir("phase1", "poisson_results.jld2"))"
