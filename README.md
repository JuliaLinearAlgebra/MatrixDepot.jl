
# ![logo](doc/logo2.png) Matrix Depot 

[![Build Status](https://travis-ci.org/weijianzhang/MatrixDepot.jl.svg?branch=master)](https://travis-ci.org/weijianzhang/MatrixDepot.jl)
| Julia 0.3 [![MatrixDepot](http://pkg.julialang.org/badges/MatrixDepot_release.svg)](http://pkg.julialang.org/?pkg=MatrixDepot&ver=release)
| Julia 0.4 [![MatrixDepot](http://pkg.julialang.org/badges/MatrixDepot_nightly.svg)](http://pkg.julialang.org/?pkg=MatrixDepot&ver=nightly)

A test matrix collection for Julia. See [Matrices](http://matrixdepotjl.readthedocs.org/en/latest/matrices.html#matrices)
for all the matrices in the collection. 

* [Documentation](http://matrixdepotjl.readthedocs.org/en/latest/)

* [Examples](http://nbviewer.ipython.org/github/weijianzhang/MatrixDepot.jl/blob/master/doc/juliadoc.ipynb)

* [Release Notes](https://github.com/weijianzhang/MatrixDepot.jl/blob/master/NEWS.md)

## Install

To install the release version, type

```julia
julia> Pkg.add("MatrixDepot")
```

## Basic Usage

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

We can type the matrix name to see the paramter options and matrix
properties.

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

From the information given by `matrixdepot("hilb")`, we observe that we
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

## References

- Nicholas J. Higham,
  "Algorithm 694, A Collection of Test Matrices in MATLAB",
  *ACM Trans. Math. Software*,
  vol. 17. (1991), pp 289-305
  [[pdf]](http://www.maths.manchester.ac.uk/~higham/narep/narep172.pdf)
  [[doi]](https://dx.doi.org/10.1145/114697.116805)
