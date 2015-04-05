
# ![logo](doc/logo2.png) Matrix Depot 

[![Build Status](https://travis-ci.org/weijianzhang/MatrixDepot.jl.svg?branch=master)](https://travis-ci.org/weijianzhang/MatrixDepot.jl)
| Julia 0.3 [![MatrixDepot](http://pkg.julialang.org/badges/MatrixDepot_release.svg)](http://pkg.julialang.org/?pkg=MatrixDepot&ver=release)
| Julia 0.4 [![MatrixDepot](http://pkg.julialang.org/badges/MatrixDepot_nightly.svg)](http://pkg.julialang.org/?pkg=MatrixDepot&ver=nightly)

An extensible test matrix collection for Julia.

* [Documentation](http://matrixdepotjl.readthedocs.org/en/latest/)

* [Demo](http://nbviewer.ipython.org/github/weijianzhang/MatrixDepot.jl/blob/master/doc/MatrixDepot_Demo.ipynb)

* [Release Notes](https://github.com/weijianzhang/MatrixDepot.jl/blob/master/NEWS.md)

**NOTE:** If you use Windows, you need to install MinGW/MSYS or
  Cygwin in order to use the UF sparse matrix collection interface.

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

We can type the matrix name to see the parameter options and matrix
properties.

```julia
julia> matrixdepot("hilb")
Hilbert matrix: 
             
 Input options: 
             
 [type,] dim: the dimension of the matrix
             
 [type,] row_dim, col_dim: the row and column dimension 
             
 ['inverse', 'ill-cond', 'symmetric', 'pos-def']
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
18-element Array{ASCIIString,1}:
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
 "lehmer"   
 "randcorr" 
 "poisson"  
 "wilkinson"
 "kms"      
 "wathen" 
```

## Interface to the UF Sparse Matrix Collection 

Use ``MatrixDepot.get(NAME)``, where ``NAME`` is ``collection
name + '/' + matrix name``, to download a test matrix from the University of
Florida Sparse Matrix Collection:
http://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html.  For
example:

```julia
julia> MatrixDepot.get("HB/1138_bus")
```
When download is complete, we can check matrix information using

```julia
julia> matrixdepot("HB/1138_bus")
%%MatrixMarket matrix coordinate real symmetric
%----------------------------------------------------------------------
% UF Sparse Matrix Collection, Tim Davis
% http://www.cise.ufl.edu/research/sparse/matrices/HB/1138_bus
% name: HB/1138_bus
% [S ADMITTANCE MATRIX 1138 BUS POWER SYSTEM, D.J.TYLAVSKY, JULY 1985.]
% id: 1
% date: 1985
% author: D. Tylavsky
% ed: I. Duff, R. Grimes, J. Lewis
% fields: title A name id date author ed kind
% kind: power network problem
%---------------------------------------------------------------------
```
and generate it with the Symbol ``:r``.

```julia
julia> matrixdepot("HB/1138_bus", :r)
1138x1138 Base.LinAlg.Symmetric{Float64,Base.SparseMatrix.SparseMatrixCSC{Float64,Int64}}:
 1474.78      0.0       0.0     …   0.0       0.0         0.0    0.0  
    0.0       9.13665   0.0         0.0       0.0         0.0    0.0  
    0.0       0.0      69.6147      0.0       0.0         0.0    0.0  
    0.0       0.0       0.0         0.0       0.0         0.0    0.0  
   -9.01713   0.0       0.0         0.0       0.0         0.0    0.0  
    0.0       0.0       0.0     …   0.0       0.0         0.0    0.0  
    0.0       0.0       0.0         0.0       0.0         0.0    0.0  
    0.0       0.0       0.0         0.0       0.0         0.0    0.0  
    0.0       0.0       0.0         0.0       0.0         0.0    0.0  
    0.0      -3.40599   0.0         0.0       0.0         0.0    0.0  
    ⋮                           ⋱             ⋮                       
    0.0       0.0       0.0         0.0       0.0         0.0    0.0  
    0.0       0.0       0.0     …   0.0     -24.3902      0.0    0.0  
    0.0       0.0       0.0         0.0       0.0         0.0    0.0  
    0.0       0.0       0.0         0.0       0.0         0.0    0.0  
    0.0       0.0       0.0         0.0       0.0         0.0    0.0  
    0.0       0.0       0.0        26.5639    0.0         0.0    0.0  
    0.0       0.0       0.0     …   0.0      46.1767      0.0    0.0  
    0.0       0.0       0.0         0.0       0.0     10000.0    0.0  
    0.0       0.0       0.0         0.0       0.0         0.0  117.647
```

The NIST Matrix Market interface is currently suspended.


See documentation for more details.

## References

- Nicholas J. Higham,
  "Algorithm 694, A Collection of Test Matrices in MATLAB",
  *ACM Trans. Math. Software*,
  vol. 17. (1991), pp 289-305
  [[pdf]](http://www.maths.manchester.ac.uk/~higham/narep/narep172.pdf)
  [[doi]](https://dx.doi.org/10.1145/114697.116805)

- T.A. Davis and Y. Hu,
  "The University of Florida Sparse Matrix Collection",
  *ACM Transaction on Mathematical Software*,
  vol. 38, Issue 1, (2011), pp 1:1-1:25
  [[pdf]](http://www.cise.ufl.edu/research/sparse/techreports/matrices.pdf)

- R.F. Boisvert, R. Pozo, K. A. Remington, R. F. Barrett, & J. Dongarra,
  " Matrix Market: a web resource for test matrix collections",
  *Quality of Numerical Software* (1996) (pp. 125-137).
  [[pdf]](ftp://ftp.idsa.prd.fr/pub/mirrors/netlib/utk/people/JackDongarra/pdf/matrixmarket.pdf)
