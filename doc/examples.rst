.. _examples:

Examples
======== 

Demo
---- 

`IJulia Notebook`_

.. _IJulia Notebook: https://github.com/weijianzhang/MatrixDepot.jl/blob/master/doc/MatrixDepot_Demo.ipynb

Getting Started
---------------

To see all the matrices in the collection, type

.. code:: 
   
   julia> matrixdepot()

   Matrices:
     1) baart            2) binomial         3) cauchy           4) chebspec      
     5) chow             6) circul           7) clement          8) deriv2        
     9) dingdong        10) fiedler         11) forsythe        12) foxgood       
    13) frank           14) grcar           15) hadamard        16) heat          
    17) hilb            18) invhilb         19) invol           20) kahan         
    21) kms             22) lehmer          23) lotkin          24) magic         
    25) minij           26) moler           27) neumann         28) oscillate     
    29) parter          30) pascal          31) pei             32) phillips      
    33) poisson         34) prolate         35) randcorr        36) rando         
    37) randsvd         38) rohess          39) rosser          40) sampling      
    41) shaw            42) toeplitz        43) tridiag         44) triw          
    45) vand            46) wathen          47) wilkinson       48) wing          

   Groups:
          data       eigen    ill-cond     inverse
       pos-def      random     regprob      sparse
     symmetric

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

    matrixdepot("hilb")
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

   matrixdepot("hadamard")
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
rectanglular Hilbert matrix by

.. code::

    matrixdepot("hilb", 4, 6)

    4x6 Array{Float64,2}:
     1.0       0.5       0.333333  0.25      0.2       0.166667
     0.5       0.333333  0.25      0.2       0.166667  0.142857
     0.333333  0.25      0.2       0.166667  0.142857  0.125   
     0.25      0.2       0.166667  0.142857  0.125     0.111111



We can aslo specify the data type

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

    matrixdepot("symmetric")

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

    matrixdepot("symmetric", "ill-cond")

    7-element Array{ASCIIString,1}:
     "hilb"   
     "cauchy" 
     "invhilb"
     "moler"  
     "pascal" 
     "pei"    
     "tridiag"



.. code:: 

    matrixdepot("inverse", "ill-cond", "symmetric")

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

We can add new groups to MatrixDepot. Since each group in Matrix
Depot is a list of strings, you can simply do, for example,

.. code:: 

    spd = matrixdepot("symmetric", "pos-def")


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


To add a group of matrices permanently for future use, we put the macro
``@addgroup`` at the beginning.

.. code:: 

    @addgroup myfav = ["lehmer", "cauchy", "hilb"]
    87

    @addgroup test_for_paper2 = ["tridiag", "sampling", "wing"]
    138

We need to reload Julia to see the changes. Type

.. code:: 

    workspace()
   
    using MatrixDepot
    matrixdepot()

    Matrices:
     1) baart            2) binomial         3) cauchy           4) chebspec      
     5) chow             6) circul           7) clement          8) deriv2        
     9) dingdong        10) fiedler         11) forsythe        12) foxgood       
    13) frank           14) grcar           15) hadamard        16) heat          
    17) hilb            18) invhilb         19) invol           20) kahan         
    21) kms             22) lehmer          23) lotkin          24) magic         
    25) minij           26) moler           27) neumann         28) oscillate     
    29) parter          30) pascal          31) pei             32) phillips      
    33) poisson         34) prolate         35) randcorr        36) rando         
    37) randsvd         38) rohess          39) rosser          40) sampling      
    41) shaw            42) toeplitz        43) tridiag         44) triw          
    45) vand            46) wathen          47) wilkinson       48) wing          

   Groups:
    data          eigen         ill-cond      inverse     
    pos-def       random        regprob       sparse      
    symmetric     myfav         test_for_paper2

    
Notice new defined groups have been included. We can use them as

.. code:: 

    matrixdepot("myfav")
    3-element Array{ASCIIString,1}:
     "lehmer"
     "cauchy"
     "hilb"  


We can remove a group using the macro ``@rmgroup``. As before, we
need to reload Julia to see the changes.

.. code:: 

    @rmproperty myfav

    153

.. code:: 

   > matrixdepot()

   Matrices:
      1) baart            2) binomial         3) cauchy           4) chebspec      
      5) chow             6) circul           7) clement          8) deriv2        
      9) dingdong        10) fiedler         11) forsythe        12) foxgood       
     13) frank           14) grcar           15) hadamard        16) heat          
     17) hilb            18) invhilb         19) invol           20) kahan         
     21) kms             22) lehmer          23) lotkin          24) magic         
     25) minij           26) moler           27) neumann         28) oscillate     
     29) parter          30) pascal          31) pei             32) phillips      
     33) poisson         34) prolate         35) randcorr        36) rando         
     37) randsvd         38) rohess          39) rosser          40) sampling      
     41) shaw            42) toeplitz        43) tridiag         44) triw          
     45) vand            46) wathen          47) wilkinson       48) wing          

  Groups:
    data          eigen         ill-cond      inverse     
    pos-def       random        regprob       sparse      
    symmetric     test_for_paper2



More Examples
-------------

An interesting test matrix is magic square. It can be generated as

.. code:: 

    M = matrixdepot("magic", 5)

    5x5 Array{Int64,2}:
     17  24   1   8  15
     23   5   7  14  16
      4   6  13  20  22
     10  12  19  21   3
     11  18  25   2   9



.. code:: 

    sum(M,1)

    1x5 Array{Int64,2}:
     65  65  65  65  65



.. code:: 

    sum(M,2)

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

    P = matrixdepot("pascal", 6)

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

    chol(P)

    6x6 Array{Float64,2}:
     1.0  1.0  1.0  1.0  1.0   1.0
     0.0  1.0  2.0  3.0  4.0   5.0
     0.0  0.0  1.0  3.0  6.0  10.0
     0.0  0.0  0.0  1.0  4.0  10.0
     0.0  0.0  0.0  0.0  1.0   5.0
     0.0  0.0  0.0  0.0  0.0   1.0


