# Formulations

## A-Formulation (Magnetostatic)

The magnetic vector potential formulation solves:

```math
\nabla \times \left(\frac{1}{\mu} \nabla \times \mathbf{A}\right) = \mathbf{J}_s
```

In 2D, ``\mathbf{A} = A_z \hat{\mathbf{z}}`` reduces to a scalar Poisson-type equation.

### Weak Form

Find ``\mathbf{A} \in V`` such that:

```math
\int_\Omega \frac{1}{\mu} (\nabla \times \mathbf{A}) \cdot (\nabla \times \mathbf{v}) \, d\Omega = \int_\Omega \mathbf{J}_s \cdot \mathbf{v} \, d\Omega \quad \forall \mathbf{v} \in V
```

## T-A Formulation (HTS)

The coupled T-A formulation solves two equations:

1. **A-equation** (full domain):
```math
\nabla \times \left(\frac{1}{\mu} \nabla \times \mathbf{A}\right) = \mathbf{J}
```

2. **T-equation** (superconductor):
```math
\nabla \cdot \mathbf{E}(\mathbf{J}) = 0
```

where ``\mathbf{J} = -\nabla T`` and the power-law relates:

```math
\mathbf{E} = E_c \left(\frac{|\mathbf{J}|}{J_c}\right)^n \frac{\mathbf{J}}{|\mathbf{J}|}
```

## Gauge Conditions

The Coulomb gauge ``\nabla \cdot \mathbf{A} = 0`` ensures uniqueness of the vector potential and is implemented via a Lagrange multiplier.
