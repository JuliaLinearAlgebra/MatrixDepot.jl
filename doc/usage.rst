Usage
-----

Every matrix in the collection is represented by a string
``"matrix_name"``, for example, the Cauchy matrix is represented by
``"cauchy"`` and the Hilbert matrix is represented by ``"hilb"``.

The matrix groups are also symbolized by strings. 
For example, the class of the symmetric
matrices is symbolized by ``"symmetric"``.

.. function:: matrixdepot()

  Return a list of all the matrices in the collection::

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
     53) wathen          54) wilkinson       55) wing          
    Groups:
     all           data          eigen         ill-cond    
     inverse       pos-def       random        regprob     
     sparse        symmetric 

.. function:: matrixdepot(matrix_name, p1, p2, ...)

  Return a matrix specified by the query string ``matrix_name``. ``p1,
  p2, ...`` are input parameters depending on ``matrix_name``. For
  example::

    julia> matrixdepot("hilb", 5, 4)
    5x4 Array{Float64,2}:
    1.0       0.5       0.333333  0.25    
    0.5       0.333333  0.25      0.2     
    0.333333  0.25      0.2       0.166667
    0.25      0.2       0.166667  0.142857
    0.2       0.166667  0.142857  0.125  

.. function:: matrixdepot(matrix_name)

  Return the documentation of ``matrix_name``, including input
  options, groups and reference. For example::

   julia> matrixdepot("moler")
      Moler Matrix
     ≡≡≡≡≡≡≡≡≡≡≡≡≡≡

   The Moler matrix is a symmetric positive definite matrix. It has one small
   eigenvalue.

   Input options:

     •  [type,] dim, alpha: dim is the dimension of the matrix, alpha is a
       scalar;

     •  [type,] dim: alpha = -1.

   Groups: ["inverse", "ill-cond", "symmetric", "pos-def"]

   References: 

   J.C. Nash, Compact Numerical Methods for Computers: Linear Algebra and
   Function Minimisation, second edition, Adam Hilger, Bristol, 1990 
   (Appendix 1).


.. function:: matrixdepot(group_name)

  Return a list of matrices which belong to group ``group_name``. For
  example::

    julia> matrixdepot("pos-def")
    11-element Array{ASCIIString,1}:
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
    "poisson"

.. function:: matrixdepot(group1, group2, ...)

  Return a list of matrices which belong to ``group1`` and ``group2``, etc. 
  For example::

    julia> matrixdepot("symmetric", "inverse", "ill-cond", "pos-def")
    7-element Array{ASCIIString,1}:
    "hilb"   
    "cauchy" 
    "invhilb"
    "moler"  
    "pascal" 
    "pei"    
    "tridiag"

.. function:: matrixdepot(num)

   Access matrix by number. For example::
     
     julia> matrixdepot(3)
     "chebspec"

.. function:: matrixdepot(num1:num2)

   Access matrix by range. For example::

     julia> matrixdepot(3:12)
     10-element Array{ASCIIString,1}:
     "chebspec"
     "chow"    
     "circul"  
     "clement" 
     "dingdong"
     "fiedler" 
     "forsythe"
     "frank"   
     "grcar"   
     "hadamard"

.. function:: matrixdepot(num, num1:num2...)

   Access matrix by a mixture of numbers and ranges. For example::

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


.. function:: matrixdepot(name, :get)

    Download a matrix from test matrix collections, where
    ``name`` is a string of collection name + ``/`` + matrix name. 
    For example::
      
      julia> matrixdepot("HB/1138_bus", :get)


.. function:: MatrixDepot.update()

    Update matrix collection database from the web server.


.. function:: matrixdepot(name)

    Output matrix information, where ``name`` is a matrix data.

.. function:: matrixdepot(name, :read)

    Generate the matrix data given by ``name``.

We can define our own groups using the macro ``@addgroup`` and
remove a defined group using ``@rmgroup``.

.. function:: @addgroup group_name = ["matrix1", "matrix2", "matrix3"]

   Create a new group ``"group_name"`` such that it has members
   ``"matrix1"``, ``"matrix2"`` and ``"matrix3"``.

.. function:: @rmgroup group_name
  
   Delete a created group ``group_name``.
