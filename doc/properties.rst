
.. _properties:

Groups
======

Groups are lists of matrix names and we use them to
categorize matrices in Matrix Depot. The list below shows
all the predefined groups in Matrix Depot and we can extend
this list by defining new groups.

Predefined Groups
-----------------

.. glossary::
   :sorted:

   symmetric
      The matrix is symmetric for some parameter values.

   inverse
      The inverse of the matrix is known explicitly.

   ill-cond
       The matrix is ill-conditioned for some parameter values.

   pos-def
       The matrix is positive definite for some parameter values.

   eigen
       Part of the eigensystem of the matrix is explicitly known.

   sparse
      The matrix is sparse.

   random
      The matrix has random entries.

   data
      The matrix has been downloaded from UF sparse matrix collection or
      the Matrix Market collection. 

   regprob
      The output is a test problem for Regularization Methods.

   all
      All the matrices in the collection. 


Adding New Groups
-----------------

New groups can be added with the macro ``@addgroup``::

    @addgroup myfav = ["lehmer", "cauchy", "hilb"]
    87

    @addgroup test_for_paper2 = ["tridiag", "sampling", "wing"]
    138

    workspace()
    
    using MatrixDepot # reload the package 
    matrixdepot()

    Matrices:
      1) baart            2) binomial         3) cauchy           4) chebspec      
      5) chow             6) circul           7)  clement          8) deriv2        
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


    matrixdepot("myfav")
    3-element Array{ASCIIString,1}:
     "lehmer"
     "cauchy"
     "hilb"
