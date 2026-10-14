# Examples

```@setup ex
using MatrixDepot, LinearAlgebra
```

## Getting started

Generate a Hilbert matrix and a circulant matrix:

```@example ex
matrixdepot("hilb", 4)
```

```@example ex
matrixdepot("circul", 5)
```

Rectangular matrices and other element types are supported where the generator allows:

```@example ex
matrixdepot("hilb", Float16, 5, 3)
```

Select matrices by group:

```@example ex
mdlist(:inverse & :illcond & :symmetric)
```

and use them in tests:

```@example ex
for name in ["lehmer", "cauchy", "hilb"]
    A = matrixdepot(name, 6)
    F = lu(A)
    println("1-norm error for $name: ", norm(A[F.p, :] - F.L * F.U, 1))
end
```

## Magic squares

```@example ex
M = matrixdepot("magic", Int, 5)
```

All row sums, column sums and both diagonal sums are equal:

```@example ex
sum(M, dims=1), sum(M, dims=2)', sum(diag(M)), sum(diag(M[:, end:-1:1]))
```

## Pascal matrices

The Cholesky factor of the Pascal matrix contains the rows of Pascal's triangle:

```@example ex
P = matrixdepot("pascal", Int, 6)
cholesky(P).U
```
