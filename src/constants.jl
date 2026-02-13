# ──────────────────────────────────────────────
# Physical Constants for HTS Electromagnetics
# ──────────────────────────────────────────────

"""Vacuum permeability μ₀ [H/m]"""
const MU_0 = 4π * 1e-7

"""Critical electric field criterion E_c [V/m] (standard for HTS characterization)"""
const E_C_DEFAULT = 1e-4

"""Boltzmann constant k_B [J/K]"""
const K_BOLTZMANN = 1.380649e-23

"""Small regularization parameter to avoid division by zero in power-law"""
const EPS_REG = 1e-20

"""Default power-law exponent for REBCO tapes"""
const N_DEFAULT = 25

"""Default critical current density J_c for REBCO [A/m²]"""
const JC_DEFAULT = 1e9
