# User Generators

Matrix Depot can be extended with your own matrix generators. `include_generator`
registers a generator and `publish_user_generators` updates the database afterwards.
`include_generator` must be called at run time, so packages defining custom generators
should call it from their `__init__` function to avoid
[precompilation issues](https://docs.julialang.org/en/v1/manual/modules/#Module-initialization-and-precompilation).
[MatrixDepotTemplate.jl](https://github.com/KlausC/MatrixDepotTemplate.jl) is a complete
example of such a package.

```@docs
MatrixDepot.include_generator
MatrixDepot.publish_user_generators
```

## Example

```@example user
using MatrixDepot, LinearAlgebra

"""
random symmetric matrix
=======================

*Input options:*

+ n: the dimension of the matrix
"""
randsym(n) = Symmetric(randn(n, n)) |> Matrix

"""
random orthogonal matrix
========================

*Input options:*

+ n: the dimension of the matrix
"""
randorth(n) = Matrix(qr(randn(n, n)).Q)

MatrixDepot.include_generator(MatrixDepot.FunctionName, "randsym", randsym)
MatrixDepot.include_generator(MatrixDepot.FunctionName, "randorth", randorth)
MatrixDepot.include_generator(MatrixDepot.Group, :random, randsym)
MatrixDepot.include_generator(MatrixDepot.Group, :symmetric, randsym)
MatrixDepot.include_generator(MatrixDepot.Group, :random, randorth)
MatrixDepot.publish_user_generators()

listnames(user(:))
```

The new generators behave like the built-in ones:

```@example user
mdinfo("randsym")
```

```@example user
A = matrixdepot("randorth", 4)
A' * A ≈ I
```

```@example user
listnames(:random)
```
