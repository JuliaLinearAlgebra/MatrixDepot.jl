# Matrix Depot

Matrix Depot is an extensible test matrix collection for Julia. It gives access to a
wealth of sample and test matrices and accompanying data:

* matrices **generated locally**, with arguments controlling the special case, and
* matrices **downloaded** from the [SuiteSparse Matrix Collection](https://sparse.tamu.edu)
  (formerly University of Florida Matrix Collection) and the
  [Matrix Market Collection](https://math.nist.gov/MatrixMarket).

Remote matrices and their metadata are cached on the local disk.

* [Source on GitHub](https://github.com/JuliaLinearAlgebra/MatrixDepot.jl)
* [Release notes](https://github.com/JuliaLinearAlgebra/MatrixDepot.jl/blob/master/NEWS.md)

## Installation

```julia-repl
pkg> add MatrixDepot
```

On Windows, MinGW/MSYS or Cygwin is needed to use the SuiteSparse and Matrix Market
interfaces.

## Quick start

```julia
using MatrixDepot

A = matrixdepot("hilb", 10)       # locally generated 10×10 Hilbert matrix
A = matrixdepot("HB/1138_bus")    # named matrix of the SuiteSparse collection
A = matrixdepot(sp(1))            # same matrix by numerical id

md = mdopen("*/bfly")             # matrix with extra data
A = md.A
co = md.coord
tx = md("Gname_10.txt")

mdinfo("gravity")                 # documentation of a matrix
md = mdopen("gravity", 10, false) # local problem with right-hand side and solution
A, b, x = md.A, md.b, md.x
```

```@contents
Pages = ["basics.md", "remote.md", "matrices.md", "regu.md", "graphs.md",
         "groups.md", "user.md", "examples.md", "api.md"]
Depth = 1
```
