# API Reference

## Main Driver

```@docs
GridapHTS.main
```

## Parameters

```@docs
GridapHTS.default_params
GridapHTS.add_default_params
GridapHTS.validate_params
```

## Materials

```@docs
GridapHTS.AbstractHTSMaterial
GridapHTS.PowerLawMaterial
GridapHTS.KimModel
GridapHTS.FieldDependentMaterial
GridapHTS.resistivity
GridapHTS.electric_field
GridapHTS.critical_current_density
GridapHTS.rebco_default
GridapHTS.bscco_default
```

## Formulations

```@docs
GridapHTS.setup_a_formulation
GridapHTS.solve_a_formulation
GridapHTS.setup_ta_formulation
GridapHTS.solve_ta_formulation
```

## Weak Forms

```@docs
GridapHTS.a_bilinear_form
GridapHTS.a_linear_form
GridapHTS.ta_residual
GridapHTS.ta_jacobian
```

## Solvers

```@docs
GridapHTS.solve_newton
GridapHTS.solve_with_continuation
```

## Gauge Conditions

```@docs
GridapHTS.apply_coulomb_gauge
```

## Mesh Utilities

```@docs
GridapHTS.load_gmsh_model
GridapHTS.setup_cartesian_model
```

## Applications

```@docs
GridapHTS.setup_iesl_benchmark
```

## Utilities

```@docs
GridapHTS.norm_safe
GridapHTS.get_nested
GridapHTS.merge_nested
GridapHTS.ensure_directory
```
