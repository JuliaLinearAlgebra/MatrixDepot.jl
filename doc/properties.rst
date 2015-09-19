
.. _properties:

Matrix Properties (Groups)
===========================


Groups
------

.. glossary::

   symmetric
      The matrix is symmetric for some parameter values.

   inverse
      The inverse of the matrix is known explicitly.

   ill-cond
       The matrix is ill-conditioned for some parameter values.

   pos-def
       The matrix is positive definite for some parameter values.

   eigen
       The eigensystem of the matrix has some known results.

   sparse
      The matrix is sparse.

   random
      The matrix has random entries.

   data
      The matrix is downloaded from UF sparse matrix collection or
      NIST Matrix Market. 

   regprob
      The output is a test problem for Regularization Methods


Examples
--------

Type ``matrixdepot("prop1", "prop2", ...)`` to see all the matrices with 
property ``prop1``, ``prop2`` etc. For example::

  julia> matrixdepot("symmetric")
  22-element Array{ASCIIString,1}:
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
  ⋮          
  "randcorr" 
  "poisson"  
  "wilkinson"
  "kms"      
  "wathen"   
  "oscillate"
  "prolate"  
  "deriv2"   
  "shaw"     

  julia> matrixdepot("symmetric", "eigen")
  11-element Array{ASCIIString,1}:
  "circul"   
  "dingdong" 
  "pascal"   
  "clement"  
  "fiedler"  
  "minij"    
  "tridiag"  
  "poisson"  
  "wilkinson"
  "wathen"   
  "oscillate"


Add New Groups
--------------

New groups can be added with the macro ``@addgroup``::

    @addgroup myfav = ["lehmer", "cauchy", "hilb"]
    87

    @addgroup test_for_paper2 = ["tridiag", "sampling", "wing"]
    138

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
