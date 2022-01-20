.. _examples:

Examples
======== 

Demo
---- 

`IJulia Notebook`_

.. _IJulia Notebook: https://github.com/JuliaMatrices/MatrixDepot.jl/blob/master/doc/MatrixDepot_Demo.ipynb

Getting Started
---------------

To see all the matrices in the collection, type

.. code:: 
   
    julia> using MatrixDepot

    julia> mdinfo()
      Currently loaded Matrices
      –––––––––––––––––––––––––––

    builtin(#)                                                                             
    ––––––––––– ––––––––––– ––––––––––– –––––––––––– ––––––––––– ––––––––––––– ––––––––––––
    1 baart     10 deriv2   19 gravity  28 kms       37 parter   46 rohess     55 ursell   
    2 binomial  11 dingdong 20 grcar    29 lehmer    38 pascal   47 rosser     56 vand     
    3 blur      12 erdrey   21 hadamard 30 lotkin    39 pei      48 sampling   57 wathen   
    4 cauchy    13 fiedler  22 hankel   31 magic     40 phillips 49 shaw       58 wilkinson
    5 chebspec  14 forsythe 23 heat     32 minij     41 poisson  50 smallworld 59 wing     
    6 chow      15 foxgood  24 hilb     33 moler     42 prolate  51 spikes                 
    7 circul    16 frank    25 invhilb  34 neumann   43 randcorr 52 toeplitz               
    8 clement   17 gilbert  26 invol    35 oscillate 44 rando    53 tridiag                
    9 companion 18 golub    27 kahan    36 parallax  45 randsvd  54 triw                   

    user(#)
    –––––––

    Groups                                                   
    ––––––– ––––– ––––– ––––––– –––––– ––––––– –––––––––     
    all     local eigen illcond posdef regprob symmetric     
    builtin user  graph inverse random sparse                

    Suite Sparse of  
    –––––––––––– ––––
    2772         2833

    MatrixMarket of 
    –––––––––––– –––
    488          498

We can generate a Hilbert matrix of size 4 by typing

.. code:: 

    matrixdepot("hilb", 4)

    4x4 Array{Float64,2}:
     1.0       0.5       0.333333  0.25    
     0.5       0.333333  0.25      0.2     
     0.333333  0.25      0.2       0.166667
     0.25      0.2       0.166667  0.142857



and generate a circul matrix of size 5 by

.. code:: 

    matrixdepot("circul", 5)

    5x5 Array{Float64,2}:
     1.0  2.0  3.0  4.0  5.0
     5.0  1.0  2.0  3.0  4.0
     4.0  5.0  1.0  2.0  3.0
     3.0  4.0  5.0  1.0  2.0
     2.0  3.0  4.0  5.0  1.0



We can type the matrix name to get help.

.. code:: 

    mdinfo("hilb")
     Hilbert matrix
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    The Hilbert matrix is a very ill conditioned matrix. It is symmetric
    positive definite and totally positive. 

    Input options:

      •  [type,] dim: the dimension of the matrix;

      •  [type,] row_dim, col_dim: the row and column dimensions.

    Groups: ["inverse", "ill-cond", "symmetric", "pos-def"]

    Reference: M. D. Choi, Tricks or treats with the Hilbert matrix, Amer. Math.
    Monthly, 90 (1983), pp. 301-312.

    N. J. Higham, Accuracy and Stability of Numerical Algorithms, Society for
    Industrial and Applied Mathematics, Philadelphia, PA, USA, 2002; sec. 28.1

.. code::

   mdinfo("hadamard")
     Hadamard matrix
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

   The Hadamard matrix is a square matrix whose entries are 1 or -1. It was
   named after Jacques Hadamard. The rows of a Hadamard matrix are orthogonal.

   Input options:

     •  [type,] n: the dimension of the matrix, n is a power of 2.

   Groups: ["inverse", "orthogonal", "eigen"]

   Reference: S. W. Golomb and L. D. Baumert, The search for Hadamard matrices,
   Amer. Math. Monthly, 70 (1963) pp. 12-17


From the information given, we can create a 4-by-6
rectangular Hilbert matrix by

.. code::

    matrixdepot("hilb", 4, 6)

    4x6 Array{Float64,2}:
     1.0       0.5       0.333333  0.25      0.2       0.166667
     0.5       0.333333  0.25      0.2       0.166667  0.142857
     0.333333  0.25      0.2       0.166667  0.142857  0.125   
     0.25      0.2       0.166667  0.142857  0.125     0.111111



We can also specify the data type

.. code:: 

    matrixdepot("hilb", Float16, 5, 3)

    5x3 Array{Float16,2}:
     1.0      0.5      0.33325
     0.5      0.33325  0.25   
     0.33325  0.25     0.19995
     0.25     0.19995  0.16663
     0.19995  0.16663  0.14282



Matrices can be accessed by groups. 

.. code:: 

    mdlist(:symmetric)

   19-element Array{ASCIIString,1}:
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
    "randsvd"  
    "kms"      
    "wathen" 

.. code:: 

    mdlist(:symmetric & :illcond)

    7-element Array{ASCIIString,1}:
     "hilb"   
     "cauchy" 
     "invhilb"
     "moler"  
     "pascal" 
     "pei"    
     "tridiag"



.. code:: 

    mdlist(:inverse & :illcond & :symmetric)

    7-element Array{ASCIIString,1}:
     "hilb"   
     "cauchy" 
     "invhilb"
     "moler"  
     "pascal" 
     "pei"    
     "tridiag"



User Defined Groups
-------------------

We can add new groups to MatrixDepot. While the predefined groups are a list of strings, you can
define user groups with with contrived patterns:

.. code:: 

    spd = mdlist(:symmetric & :posdef) # intersection of two groups - also `|` and `~` are supported, see ?MatrixDepot


    10-element Array{ASCIIString,1}:
     "hilb"   
     "cauchy" 
     "circul" 
     "invhilb"
     "moler"  
     "pascal" 
     "pei"    
     "minij"  
     "tridiag"
     "lehmer" 



.. code::

    myprop = ["lehmer", "cauchy", "hilb"]

    3-element Array{ASCIIString,1}:
     "lehmer"
     "cauchy"
     "hilb"  



Then use it in your tests like

.. code:: 

    for matrix in myprop
        A = matrixdepot(matrix, 6)
        L, U, p = lu(A) #LU factorization
        err = norm(A[p,:] - L*U, 1) # 1-norm error
        println("1-norm error for $matrix matrix is ", err)
    end    

    1-norm error for lehmer matrix is 1.1102230246251565e-16
    1-norm error for cauchy matrix is 5.551115123125783e-17
    1-norm error for hilb matrix is 2.7755575615628914e-17


To add a group of matrices permanently for future use, we call
`setgroup!` at the beginning. Note the `:` in the group name.

.. code:: 

    setgroup!(:myfav, ["lehmer", "cauchy", "hilb"])

    setgroup!(:test_for_paper2, ["tridiag", "sampling", "wing"])

You can see the changes immediately:

.. code:: 

    mdinfo()
      Currently loaded Matrices
      –––––––––––––––––––––––––––

    builtin(#)                                                                             
    ––––––––––– ––––––––––– ––––––––––– –––––––––––– ––––––––––– ––––––––––––– ––––––––––––
    1 baart     10 deriv2   19 gravity  28 kms       37 parter   46 rohess     55 ursell   
    2 binomial  11 dingdong 20 grcar    29 lehmer    38 pascal   47 rosser     56 vand     
    3 blur      12 erdrey   21 hadamard 30 lotkin    39 pei      48 sampling   57 wathen   
    4 cauchy    13 fiedler  22 hankel   31 magic     40 phillips 49 shaw       58 wilkinson
    5 chebspec  14 forsythe 23 heat     32 minij     41 poisson  50 smallworld 59 wing     
    6 chow      15 foxgood  24 hilb     33 moler     42 prolate  51 spikes                 
    7 circul    16 frank    25 invhilb  34 neumann   43 randcorr 52 toeplitz               
    8 clement   17 gilbert  26 invol    35 oscillate 44 rando    53 tridiag                
    9 companion 18 golub    27 kahan    36 parallax  45 randsvd  54 triw                   

    user(#)
    –––––––

    Groups                                                                    
    ––––––– ––––– ––––– ––––––– –––––– ––––––– ––––––––– –––––––––––––––      
    all     local eigen illcond posdef regprob symmetric test_for_paper2      
    builtin user  graph inverse random sparse  myfav                          

    Suite Sparse of  
    –––––––––––– ––––
    2772         2833

    MatrixMarket of 
    –––––––––––– –––
    488          498
    
Notice new defined groups have been included. We can use them as

.. code:: 

    mdlist(:myfav)
    3-element Array{ASCIIString,1}:
     "lehmer"
     "cauchy"
     "hilb"  


We can remove a group using the macro ``@rmgroup``. As before, we
need to reload Julia to see the changes.

.. code:: 

    @rmgroup myfav

.. code:: 

   listgroups()
    14-element Array{Symbol,1}:
     :all            
     :builtin        
     :local          
     :user           
     :eigen          
     :graph          
     :illcond        
     :inverse        
     :posdef         
     :random         
     :regprob        
     :sparse         
     :symmetric      

More Examples
-------------

An interesting test matrix is magic square. It can be generated as

.. code:: 

    M = matrixdepot("magic", Int, 5)

    5x5 Array{Int64,2}:
     17  24   1   8  15
     23   5   7  14  16
      4   6  13  20  22
     10  12  19  21   3
     11  18  25   2   9



.. code:: 

    sum(M, dims=1)

    1x5 Array{Int64,2}:
     65  65  65  65  65



.. code:: 

    sum(M, dims=2)

    5x1 Array{Int64,2}:
     65
     65
     65
     65
     65



.. code:: 

    sum(diag(M))

    65



.. code:: 

    p = [5:-1:1]
    sum(diag(M[:,p]))

    65



Pascal Matrix can be generated as

.. code:: 

    P = matrixdepot("pascal", Int, 6)

    6x6 Array{Int64,2}:
     1  1   1   1    1    1
     1  2   3   4    5    6
     1  3   6  10   15   21
     1  4  10  20   35   56
     1  5  15  35   70  126
     1  6  21  56  126  252



Notice the Cholesky factor of the Pascal matrix has Pascal's triangle
rows.

.. code:: 

    cholesky(P)

    6x6 UpperTriangular{Float64,Array{Float64,2}}:
     1.0  1.0  1.0  1.0  1.0   1.0
     0.0  1.0  2.0  3.0  4.0   5.0
     0.0  0.0  1.0  3.0  6.0  10.0
     0.0  0.0  0.0  1.0  4.0  10.0
     0.0  0.0  0.0  0.0  1.0   5.0
     0.0  0.0  0.0  0.0  0.0   1.0

