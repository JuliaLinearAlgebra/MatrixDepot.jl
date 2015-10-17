
# ![logo](doc/logo2.png) Matrix Depot 

[![Build Status](https://travis-ci.org/weijianzhang/MatrixDepot.jl.svg?branch=master)](https://travis-ci.org/weijianzhang/MatrixDepot.jl)

An extensible test matrix collection for Julia.

* [Documentation](http://matrixdepotjl.readthedocs.org/en/latest/)

* [Demo](https://github.com/weijianzhang/MatrixDepot.jl/blob/master/doc/MatrixDepot_Demo.ipynb)

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

Matrices:
   1) baart            2) binomial         3) blur             4) cauchy        
   5) chebspec         6) chow             7) circul           8) clement       
   9) deriv2          10) dingdong        11) fiedler         12) forsythe      
  13) foxgood         14) frank           15) gravity         16) grcar         
  17) hadamard        18) hankel          19) heat            20) hilb          
  21) invhilb         22) invol           23) kahan           24) kms           
  25) lehmer          26) lotkin          27) magic           28) minij         
  29) moler           30) neumann         31) oscillate       32) parter        
  33) pascal          34) pei             35) phillips        36) poisson       
  37) prolate         38) randcorr        39) rando           40) randsvd       
  41) rohess          42) rosser          43) sampling        44) shaw          
  45) toeplitz        46) tridiag         47) triw            48) vand          
  49) wathen          50) wilkinson       51) wing          
Groups:
  all           data          eigen         ill-cond    
  inverse       pos-def       random        regprob     
  sparse        symmetric  
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
             
 1. [type,] dim: the dimension of the matrix;
             
 2. [type,] row_dim, col_dim: the row and column dimension.
             
 ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             
 Reference: M. D. Choi, Tricks or treats with the Hilbert matrix,
             Amer. Math. Monthly, 90 (1983), pp. 301-312.
             N. J. Higham, Accuracy and Stability of Numerical Algorithms,
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 2002; sec. 28.1.
```

We can also specify the data type

```julia
julia> matrixdepot("hilb", Float16, 5, 3)
5x3 Array{Float16,2}:
 1.0      0.5      0.33325
 0.5      0.33325  0.25   
 0.33325  0.25     0.19995
 0.25     0.19995  0.16663
 0.19995  0.16663  0.14282

julia> matrixdepot("hilb", Rational, 4)
4x4 Array{Rational{T<:Integer},2}:
 1//1  1//2  1//3  1//4
 1//2  1//3  1//4  1//5
 1//3  1//4  1//5  1//6
 1//4  1//5  1//6  1//7
```

Matrices can be accessed by number, range or a mixture of number and range.

```julia
julia> matrixdepot(5)
"chow"

julia> matrixdepot(5:10)
6-element Array{AbstractString,1}:
 "chow"    
 "circul"  
 "clement" 
 "deriv2"  
 "dingdong"
 "fiedler"

julia> matrixdepot(1:4, 6, 10:15)
11-element Array{AbstractString,1}:
 "baart"   
 "binomial"
 "cauchy"  
 "chebspec"
 "circul"  
 "fiedler" 
 "forsythe"
 "foxgood" 
 "frank"   
 "gravity" 
 "grcar" 
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

## Extend Matrix Depot

We can add more matrices to Matrix Depot by downloading them from UF
sparse matrix collection and Matrix Market. See
[here](http://matrixdepotjl.readthedocs.org/en/latest/interface.html)
for more details.

We can add new groups of matrices. See
[here](http://matrixdepotjl.readthedocs.org/en/latest/properties.html)
for more details.

## Interface to the UF Sparse Matrix Collection 

Use ``matrixdepot(NAME, :get)``, where ``NAME`` is ``collection
name + '/' + matrix name``, to download a test matrix from the University of
Florida Sparse Matrix Collection:
http://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html.  For
example:

```julia
julia> matrixdepot("HB/1138_bus", :get)
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
and generate it with the Symbol ``:r`` or ``:read``.

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

The NIST Matrix Market interface is similar. See
[documentation](http://matrixdepotjl.readthedocs.org/en/latest/interface.html#interface-to-nist-matrix-market)
for more details.


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
  [[pdf]](http://www.netlib.org/utk/people/JackDongarra/pdf/matrixmarket.pdf)

- Per Christian Hansen,
  "Test Matrices for Regularization Methods",
  *SIAM Journal on Scientific Computing*,
  vol. 16, 2, (1995) pp.506-512.
  [[pdf]](http://epubs.siam.org/doi/abs/10.1137/0916032)
  [[doi]](https://dx.doi.org/10.1137/0916032)
