
.. _properties:

Matrix Properties
=================

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


Main Properties
---------------

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

Other Properties
----------------

.. glossary::

   sparse
      The matrix is sparse.

   random
      The matrix has random entries.

   data
      Matrix data downloaded from UF sparse matrix collection or
      NIST Matrix Market. 

   regu
      Test problems for Regularization Methods
