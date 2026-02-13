# ──────────────────────────────────────────────
# Public API Exports
# ──────────────────────────────────────────────

# Main driver
export main

# Parameters
export default_params, add_default_params, validate_params

# Constants
export MU_0, E_C_DEFAULT, N_DEFAULT, JC_DEFAULT

# Materials
export AbstractHTSMaterial, PowerLawMaterial
export resistivity, electric_field, critical_current_density
export KimModel, FieldDependentMaterial
export rebco_default, bscco_default

# Formulations
export setup_a_formulation, solve_a_formulation
export setup_ta_formulation, solve_ta_formulation
export a_bilinear_form, a_linear_form
export ta_residual, ta_jacobian

# Solvers
export solve_newton, solve_with_continuation

# Gauges
export apply_coulomb_gauge

# Mesh utilities
export load_gmsh_model, setup_cartesian_model

# Boundary Conditions
export apply_dirichlet_bc, apply_neumann_bc

# Applications
export setup_iesl_benchmark

# Utilities
export norm_safe, get_nested, merge_nested, ensure_directory
