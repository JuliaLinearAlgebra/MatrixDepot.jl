# Test Problems for Regularization Methods

```@setup regu
using MatrixDepot
```

A Fredholm integral equation of the first kind (in one dimension) can be written as

```math
\int_0^1 K(s,t) f(t)\, dt = g(s), \quad 0 \leq s \leq 1,
```

where ``g`` and the kernel ``K`` are known functions and ``f`` is the unknown solution.
This is a classical linear ill-posed problem: an arbitrarily small perturbation of the
data can cause an arbitrarily large perturbation of the solution. In computerized
tomography, for example, ``K`` is an X-ray source, ``f`` is the scanned object and ``g``
is the measured damping of the X-rays.

Discretization (by quadrature or the Galerkin method) yields a linear system ``Ax = b``.
All problems of the group `:regprob` are derived this way. By default only the matrix is
generated; with `matrixonly = false` as last argument, [`mdopen`](@ref) returns a
descriptor with the matrix `A`, right-hand side `b` and solution `x`:

```@example regu
A = matrixdepot("deriv2", 4)
```

```@example regu
md = mdopen("deriv2", 4, false)
md.x
```

## Problems

### baart

```@docs
MatrixDepot.baart
```

### blur

```@docs
MatrixDepot.blur
```

### deriv2

```@docs
MatrixDepot.deriv2
```

### foxgood

```@docs
MatrixDepot.foxgood
```

### gravity

```@docs
MatrixDepot.gravity
```

### heat

```@docs
MatrixDepot.heat
```

### parallax

```@docs
MatrixDepot.parallax
```

### phillips

```@docs
MatrixDepot.phillips
```

### shaw

```@docs
MatrixDepot.shaw
```

### spikes

```@docs
MatrixDepot.spikes
```

### ursell

```@docs
MatrixDepot.ursell
```

### wing

```@docs
MatrixDepot.wing
```

