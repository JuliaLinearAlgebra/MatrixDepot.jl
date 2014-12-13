# Matrix Depot 

A test matrix collection for Julia.

[![Build Status](https://travis-ci.org/weijianzhang/MatrixDepot.jl.svg?branch=master)](https://travis-ci.org/weijianzhang/MatrixDepot.jl)

  
## Install

To install the release version, type

```julia
julia> Pkg.add("MatrixDepot")
```
To install the latest development version, type

```julia
julia> Pkg.clone("MatrixDepot")
```

## Usage

First load the package:

```julia
julia> using MatrixDepot
```

The only function will be exported is `matrixdepot`.

```julia
julia> ? matrixdepot
INFO: Loading help data...
matrixdepot (generic function with 15 methods)
```

Every matrix in the collection is represented by a string `matrix_name`, for
example, the Cauchy matrix is represented by `"cauchy"` and the Hilbert matrix
is represented by `"hilb"`.

The properties of the matrices in the collection are also symbolized by strings
`propertry_name`. For example, the class of the symmetric matrices is symbolized
by `"symmetric"`.

* `matrixdepot()` returns a list of all the matrices in the collection.

* `matrixdepot(matrix_name, p1, p2, ...)` returns a matrix specified by the
query `matrix_name`. `matrix_name` is a query string. `p1, p2, ...` are input
parameters depending on `matrix_name`.

* `matrixdepot(matrix_name)` returns the parameter options and the properties of
`matrix_name`.

* `matrixdepot(prop1, prop2, ...)` returns a list of matrices with the property
`prop1`, `prop2`, etc.

### Matrices in the Collection

| Matrices                  | Strings   | Matrices                 | Strings   |
|:-------------------------:|:---------:|:------------------------:|:---------:|
| Circul Matrix             | "circul"  | Frank Matrix             | "frank"   | 
| Hilbert Matrix            | "hilb"    | Grcar Matrix             | "grcar"   |
| Inverse of Hilbert Matrix | "invhilb" | Dingdong Matrix          | "dingdong"|
| Hadamard Matrix           | "hadamard"| Magic Matrix             | "magic"   |
| Cauchy Matrix             | "cauchy"  | Forsythe Matrix          | "forsythe"|
| Triw Matrix               | "triw"    | Moler Matrix             | "moler"   |
| Pascal Matrix             | "pascal"  | Kahan Matrix             | "kahan"   |
| Pei Matrix                | "pei"     | Vandermonde Matrix       | "vand"    |
| Involutory Matrix         | "invol"   | Cheb. spec. diff. Matrix | "chebspec"| 
| Lotkin Matrix             | "lotkin"  | Clement Matrix           | "clement" |
| Fiedler Matrix            | "fiedler" | MIN[I,J] Matrix          | "minij"   |
| Binomial Matrix           | "binomial"| Tridiagonal Matrix       | "tridiag" |

### Matrix properties in the Collection

* `"symmetric"`: the matrix is symmetric for some parameter values.

* `"inverse"`: the inverse of the matrix is known explicitly.

* `"ill-cond"`: the matrix is ill-conditioned for some parameter values.

* `"pos-def"`: the matrix is symmetric postive definite for some parameter values.

* `"eigen"`: the eigensystem of the matrix has some known results (explicit 
formulas for eigenvalues, eigenvectors, bounds of eigenvalues, etc).

## Examples

To see all the matrices in the collection, type

```julia
julia> matrixdepot()
```

We can generate a Hilbert matrix of size 4 by typing

```julia
julia> matrixdepot("hilb", 4)
4x4 Array{Float64,2}:
 1.0       0.5       0.333333  0.25    
 0.5       0.333333  0.25      0.2     
 0.333333  0.25      0.2       0.166667
 0.25      0.2       0.166667  0.142857
```

and generate a circul matrix of size 5 by

```julia
julia> matrixdepot("circul", 5)
5x5 Array{Float64,2}:
 1.0  2.0  3.0  4.0  5.0
 5.0  1.0  2.0  3.0  4.0
 4.0  5.0  1.0  2.0  3.0
 3.0  4.0  5.0  1.0  2.0
 2.0  3.0  4.0  5.0  1.0
```

If we are not sure about the paramters or properties of a matrix, we can type
the matrix name

```julia
julia> matrixdepot("hilb")
Hilbert matrix: 
                  
Input options: 
                  
(type), dim: the dimension of the matrix
                  
(type), row_dim, col_dim: the row and column dimension 
                  
['inverse', 'ill-cond', 'symmetric', 'pos-def']


julia> matrixdepot("hadamard")
Hadamard matrix: 
                  
Input options: 
                  
(type), dim: the dimension of the matrix, n is a power of 2 
                  
['inverse', 'orthogonal', 'eigen']
```

From the information given by `matrixdepot("hilb")`, we notice that we
can create a 4-by-6 rectanglular Hilbert matrix by

```julia
julia> matrixdepot("hilb", 4, 6)
4x6 Array{Float64,2}:
 1.0       0.5       0.333333  0.25      0.2       0.166667
 0.5       0.333333  0.25      0.2       0.166667  0.142857
 0.333333  0.25      0.2       0.166667  0.142857  0.125   
 0.25      0.2       0.166667  0.142857  0.125     0.111111
```

We can aslo specify the data type

```julia
julia> matrixdepot("hilb", Float16, 5, 3)
5x3 Array{Float16,2}:
 1.0      0.5      0.33325
 0.5      0.33325  0.25   
 0.33325  0.25     0.19995
 0.25     0.19995  0.16663
 0.19995  0.16663  0.14282
```

By typing a matrix name, we can see what properties that matrix have.
Conversely, if we type a property (or properties), we can see all the 
matrices (in the collection) having that property (or properties).

```julia
julia> matrixdepot("symmetric")
12-element Array{ASCIIString,1}:
 "hilb"    
 "cauchy"  
 "circul"  
 "dingdong"
 "invhilb" 
 "moler"   
 "pascal"  
 "pei"     
 "clement" 
 "fiedler" 
 "minij"   
 "tridiag" 

julia> matrixdepot("inverse", "ill-cond", "symmetric")
7-element Array{ASCIIString,1}:
 "hilb"   
 "cauchy" 
 "invhilb"
 "moler"  
 "pascal" 
 "pei"    
 "tridiag"
```  

Given a property (or properites), we can loop through all the matrices 
having this propery (or properties)

```julia
# Multiply all matrices of the class "symmetric", "ill-cond" and "inverse".
julia> A = eye(4)
julia> print("Identity matrix")
julia> for mat in matrixdepot("symmetric", "ill-cond", "inverse")
           print(" x $mat matrix")
           A = A * full(matrixdepot(mat, 4))
       end
julia> println(" =")
julia> A
   
Identity matrix x hilb matrix x cauchy matrix x invhilb matrix x moler matrix x pascal matrix x pei matrix x tridiag matrix =

4x4 Array{Float64,2}:
 153.12    -11.919    -15.4345   296.937
 109.896    -8.91857  -11.5976   214.433
  86.7524   -7.15714   -9.32857  169.702
  71.9139   -5.98707   -7.81497  140.876 
```

## Documentation

More details can be found [here](http://nbviewer.ipython.org/github/weijianzhang/MatrixDepot.jl/blob/master/doc/juliadoc.ipynb).

## References

- Nicholas J. Higham,
  "Algorithm 694, A Collection of Test Matrices in MATLAB",
  *ACM Trans. Math. Software*,
  vol. 17. (1991), pp 289-305
  [[pdf]](http://www.maths.manchester.ac.uk/~higham/narep/narep172.pdf)
  [[doi]](https://dx.doi.org/10.1145/114697.116805)