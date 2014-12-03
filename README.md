# Matrix Depot for Julia

Matrix Depot is a test matrix collection.

[![Build Status](https://travis-ci.org/weijianzhang/MatrixDepot.jl.svg?branch=master)](https://travis-ci.org/weijianzhang/MatrixDepot.jl)

  
## Install

To install the package, type

```julia
julia> Pkg.update()
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
matrixdepot (generic function with 10 methods)
```

Every matrix in the collection is represented by a string `matrix_name`, for
example, the Cauchy matrix is represented by `"cauchy"` and the Hilbert matrix
is represented by `"hilb"`.

The properties of the matrices in the collection are also symbolized by strings
`propertry_name`. For example, the class of the symmetric matrices is symbolized
by `"symmetric"`.

* `matrixdepot(matrix_name, p1, p2, ...)` returns a matrix specified by the
query `matrix_name`. `matrix_name` is a query string. `p1, p2, ...` are input
parameters depending on `matrix_name`.

* `matrixdepot(matrix_name)` returns the parameter options and the properties of
`matrix_name`.

* `matrixdepot(property_name)` returns a list of matrices with the property
`property_name`.


| Matrices                  | String    | Matrices       | String    |
|---------------------------|-----------|----------------|-----------|
| Circul Matrix             | "circul"  | Frank Matrix   | "frank"   | 
| Hilbert Matrix            | "hilb"    | Grcar Matrix   | "grcar"   |
| Inverse of Hilbert Matrix | "invhilb" | Dingdong Matrix| "dingdong"|
| Hadamard Matrix           | "hadamard"| Magic Matrix   | "magic"   |
| Cauchy Matrix             | "cauchy"  | Forsythe Matrix| "forsythe"|


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
Conversely, if we type a property, we can see all the matrices (in the
collection) having that property.

```julia
julia> matrixdepot("symmetric")
5-element Array{ASCIIString,1}:
 "hilb"  
 "cauchy"
 "circul"
 "dingdong"
 "invhilb"
```  

Given a property, we can loop through all the matrices having this propery

```julia
 # Multiply all matrices of the class "inverse"
julia> A = eye(4);
julia> print("Identity matrix")
julia> for mat in matrixdepot("inverse")
         print(" x $mat matrix")
         A *=  matrixdepot(mat, 4)    
       end
julia> println(" =")
julia> A    
Identity matrix x hilb matrix x hadamard matrix x cauchy matrix x invhilb matrix =
4x4 Array{Float64,2}:
 -0.00595238   2.20238   0.047619   1.75    
 -0.00166667   1.31667   0.0666667  0.616667
 -0.000714286  0.964286  0.052381   0.316667
 -0.00037415   0.767007  0.0401361  0.192857
```

## Documentation

More details can be found [here](http://nbviewer.ipython.org/github/weijianzhang/MatrixDepot.jl/blob/master/doc/juliadoc.ipynb).
