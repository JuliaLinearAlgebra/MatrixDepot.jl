
# ![logo](doc/logo2.png) Matrix Depot

[![Build Status](https://travis-ci.org/JuliaMatrices/MatrixDepot.jl.svg?branch=master)](https://travis-ci.org/JuliaMatrices/MatrixDepot.jl)
[![codecov.io](https://codecov.io/github/JuliaMatrices/MatrixDepot.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaMatrices/MatrixDepot.jl?branch=master)

An extensible test matrix collection for Julia.

* [Documentation](http://matrixdepotjl.readthedocs.org/en/latest/)

* [Release Notes](https://github.com/weijianzhang/MatrixDepot.jl/blob/master/NEWS.md)

**NOTE:** If you use Windows, you need to install MinGW/MSYS or
  Cygwin in order to use the UF sparse matrix collection interface.

## Install

To install the release version, type

```julia
julia> Pkg.add("MatrixDepot")
```

## New API

This is a re-factoring of the existing code.
While the vintage API `matrixdepot(...)` is still supported it adds quite some new features.
1. `list(pattern)` delivers a list of full names matching the pattern. several variants of pattern are supported:
1.1 strings with or without shell pattern characters (`'*'`, `'?'`). Double star also matches `'/'`
1.2 regular expressions
1.3 integer matrix identifiers (as provided by the UFL collection) and ranges thereof
1.4 intersections of patterns (represented as a tuple of patterns)
1.5 unions of patterns (represented by an array of patterns)
1.6 group names - given as symbols, e.g. `:illcond` instead of `"ill-cond"`.
2. `info(pattern)` delivers information about all matrices matched by pattern as markdown object
3. `matrix(pattern, args...)` delivers a generated (with arguments) or loaded problem matrix
4.  `rhs(pattern)` delivers the b-side of the problem - if provided by the remote matrix.
5. `solution(pattern)` delivers the x-side of the problem - if provided as metadata.
6. `mreader(pattern, name)` delivers other metadata.
7. `load(pattern)` loads all remote matrices matching the pattern from the dedicated servers.

The Matrixmarket-format reader has been completely re-written. It is now more reliable and faster than the previous approach, which used CHOLMOD.Sparse.

All retrieved data and metadata may be cached in memory on demand.
Information about the matrixmarket-format and the metadata are accessible using `data = mdopen(pattern)` or `data = MatrixDepot.mdata(pattern)`.
`mdclose(data)` removes all cached data.

## Basic Usage

To see an overview of the matrices in the collection, type

```julia
julia> using MatrixDepot
julia> mdinfo()
  Currently loaded Matrices
  –––––––––––––––––––––––––––

builtin(#)                                                                 
––––––––––– ––––––––––– ––––––––––– –––––––––––– ––––––––––––– ––––––––––––
1 baart     11 dingdong 21 hadamard 31 magic     41 poisson    51 spikes   
2 binomial  12 erdrey   22 hankel   32 minij     42 prolate    52 toeplitz 
3 blur      13 fiedler  23 heat     33 moler     43 randcorr   53 tridiag  
4 cauchy    14 forsythe 24 hilb     34 neumann   44 rando      54 triw     
5 chebspec  15 foxgood  25 invhilb  35 oscillate 45 randsvd    55 ursell   
6 chow      16 frank    26 invol    36 parallax  46 rohess     56 vand     
7 circul    17 gilbert  27 kahan    37 parter    47 rosser     57 wathen   
8 clement   18 golub    28 kms      38 pascal    48 sampling   58 wilkinson
9 companion 19 gravity  29 lehmer   39 pei       49 shaw       59 wing     
10 deriv2   20 grcar    30 lotkin   40 phillips  50 smallworld             

user(#)  
–––––––––
1 randsym

Groups                                                   
––––––– ––––– ––––– ––––––– –––––– ––––––– –––––––––     
all     local eigen illcond posdef regprob symmetric     
builtin user  graph inverse random sparse                

UFL / TAMU of  
–––––––––– ––––
1          2757

MatrixMarket of 
–––––––––––– –––
0            498

```

We can generate a 4-by-4 Hilbert matrix by typing

```julia
julia> matrixdepot("hilb", 4)
4x4 Array{Float64,2}:
 1.0       0.5       0.333333  0.25    
 0.5       0.333333  0.25      0.2     
 0.333333  0.25      0.2       0.166667
 0.25      0.2       0.166667  0.142857
```

We can type the matrix name to get documentation about the matrix.

```julia
julia> mdinfo("hilb")
     Hilbert matrix
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

  The Hilbert matrix has (i,j) element 1/(i+j-1). It is notorious for being
  ill conditioned. It is symmetric positive definite and totally positive.

  Input options:

    •  [type,] dim: the dimension of the matrix;

    •  [type,] row_dim, col_dim: the row and column dimensions.

  Groups: ["inverse", "ill-cond", "symmetric", "pos-def"]

  References:

  M. D. Choi, Tricks or treats with the Hilbert matrix, Amer. Math. Monthly,
  90 (1983), pp. 301-312.

  N. J. Higham, Accuracy and Stability of Numerical Algorithms, second
  edition, Society for Industrial and Applied Mathematics, Philadelphia, PA,
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

julia> matrixdepot("hilb", Rational{Int}, 4)
4x4 Array{Rational{T<:Integer},2}:
 1//1  1//2  1//3  1//4
 1//2  1//3  1//4  1//5
 1//3  1//4  1//5  1//6
 1//4  1//5  1//6  1//7
```

Matrices can be accessed by a variety of patterns and composed patterns.
Integer numbers refer to the ident numbers of the TAMU/UFl collection.

```julia
julia> mdlist(uf(1))    # here uf(1) is the ident number of the UFL collection
list(1)
–––––––––––
HB/1138_bu

julia> mdlist(builtin(1, 5:10))    # the internal numbering of the builtin-functions
list(7)
––––––– –––––––– –––– –––––– ––––––– ––––––––– ––––––
baart   chebspec chow circul clement companion deriv2

julia> MatrixDepot.list(builtin(1:4, 6, 10:15) | user(1:10) )
12-element Array{String,1}:
 "baart"
 "binomial"
 "blur"
 "cauchy"
 "chow"
 "deriv2"
 "dingdong"
 "erdrey"
 "fiedler"
 "forsythe"
 "foxgood"
 "rand
```

While the `mdlist` command renders the output as markdown table, the internal
`MatrixDepot.list` produces an array of valid matrix names.

We can type a group name to see all the matrices in that group. Group names are
always written as symbols to distinguish them form matrix names and pattern, which
are always strings.

```julia
julia> mdlist(:symmetric)
list(22)
–––––––– –––––––– ––––––– –––––– ––––––––– –––––––– ––––––– –––––––––
cauchy   dingdong hilb    lehmer oscillate poisson  randsym wilkinson
circul   fiedler  invhilb minij  pascal    prolate  tridiag
clement  hankel   kms     moler  pei       randcorr wathen



```

## Extend Matrix Depot

We can add more matrices to Matrix Depot by downloading them from UF
sparse matrix collection and Matrix Market. See
[here](http://matrixdepotjl.readthedocs.org/en/latest/interface.html)
for more details.
In addition,
we can add [new matrix generators](http://matrixdepotjl.readthedocs.org/en/latest/user.html)
and define [new groups of matrices](http://matrixdepotjl.readthedocs.org/en/latest/properties.html).


## Interface to the UF Sparse Matrix Collection

Use ``load(NAME)``, where ``NAME`` is ``collection_name + '/' + matrix_name``, to download a test matrix from the University of
Florida Sparse Matrix Collection:
http://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html.  For
example:

```julia
julia> matrixdepot("HB/1138_bus", :get)
```

Use ``matrixdepot(collection_name/*, :get )`` to download a group of all the matrices in ``collection_name`` from UF sparse matrix collection. For example:

```julia
julia> matrixdepot("MathWorks/*", :get)
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
1138x1138 Symmetric{Float64,SparseMatrixCSC{Float64,Int64}}:
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

Note ``matrixdepot()`` displays all the matrices in the collection,
including the newly downloaded matrices.

```julia
julia> matrixdepot()

Matrices:
   1) baart            2) binomial         3) blur             4) cauchy        
   5) chebspec         6) chow             7) circul           8) clement       
   9) companion       10) deriv2          11) dingdong        12) fiedler       
  13) forsythe        14) foxgood         15) frank           16) golub         
  17) gravity         18) grcar           19) hadamard        20) hankel        
  21) heat            22) hilb            23) invhilb         24) invol         
  25) kahan           26) kms             27) lehmer          28) lotkin        
  29) magic           30) minij           31) moler           32) neumann       
  33) oscillate       34) parter          35) pascal          36) pei           
  37) phillips        38) poisson         39) prolate         40) randcorr      
  41) rando           42) randsvd         43) rohess          44) rosser        
  45) sampling        46) shaw            47) spikes          48) toeplitz      
  49) tridiag         50) triw            51) ursell          52) vand          
  53) wathen          54) wilkinson       55) wing            56) HB/1138_bus   

Groups:
  all           data          eigen         ill-cond    
  inverse       pos-def       random        regprob     
  sparse        symmetric  
```

The NIST Matrix Market interface is similar. See
[documentation](http://matrixdepotjl.readthedocs.org/en/latest/interface.html#interface-to-nist-matrix-market)
for more details.


## References

- Weijian Zhang and Nicholas J. Higham,
  "Matrix Depot: An Extensible Test Matrix Collection for Julia",
  *PeerJ Comput. Sci.*, 2:e58 (2016),
  [[pdf]](https://peerj.com/articles/cs-58/)

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
