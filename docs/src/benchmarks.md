# Benchmarks

## IESL 2D Rectangular Tape

The primary validation benchmark from [htsmodelling.com](https://htsmodelling.com/).

### Problem Description

A single HTS tape cross-section (width 12 mm, thickness 1 μm) under an applied perpendicular magnetic field. The screening current distribution is computed and compared against reference solutions.

### Parameters

| Parameter | Value | Unit |
|:----------|:------|:-----|
| Tape width | 12 | mm |
| Tape thickness | 1 | μm |
| ``J_c`` | ``3 \times 10^{10}`` | A/m² |
| ``n`` | 21 | - |
| ``E_c`` | ``10^{-4}`` | V/m |
| ``B_{applied}`` | 0 - 1 | T |

### Usage

```julia
using GridapHTS
params = setup_iesl_benchmark(; jc=3e10, n=21, b_max=0.5)
xh, fullparams, info = GridapHTS.main(params)
```

## Validation Status

| Benchmark | Status | Phase |
|:----------|:-------|:------|
| Poisson equation | Implemented | 1 |
| Magnetostatic A-formulation | Implemented | 1 |
| IESL 2D Tape | Planned | 2 |
| 3D Solenoid | Planned | 3 |
