Usage
-----

Every matrix in the collection is represented by a string
``"matrix_name"``, for example, the Cauchy matrix is represented by
``"cauchy"`` and the Hilbert matrix is represented by ``"hilb"``.

The matrix groups are noted as symbols. 
For example, the class of the symmetric
matrices is symbolized by ``:symmetric``.

.. function:: mdinfo()

  Return a list of all the matrices in the collection::

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
    2777         2893

    MatrixMarket of 
    –––––––––––– –––
    488          498

.. function:: matrixdepot(matrix_name, p1, p2, ...)

  Return a matrix specified by the query string ``matrix_name``.
  The string must be a matrix name or a query pattern which matches exactly one matrix.
  ``p1, p2, ...`` are input parameters depending on ``matrix_name``.
  For example::

    julia> matrixdepot("hilb", 5, 4)
    5x4 Array{Float64,2}:
    1.0       0.5       0.333333  0.25    
    0.5       0.333333  0.25      0.2     
    0.333333  0.25      0.2       0.166667
    0.25      0.2       0.166667  0.142857
    0.2       0.166667  0.142857  0.125  

.. function:: mdinfo(matrix_name)

  Return the documentation of ``matrix_name``, including input options,
  groups and reference. For example::

   julia> mdinfo("moler")
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


.. function:: mdlist(group_name)

  Return a list of matrices which belong to group ``group_name`` (symbol) as an array.
  For example::

    julia> mdlist(:posdef)
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

.. function:: listnames(group1 & group2 & ...)

  Return a iformatted list of matrices which belong to ``group1`` and ``group2``, etc. 
  For example::

  julia> listnames(:symmetric & :inverse & :illcond & :posdef)
  list(8)
  ––––––– –––– ––––––– ––– ––––– –––––– ––– –––––––
  cauchy  hilb invhilb kms moler pascal pei tridiag

.. function:: mdlist({builtin\user\sp\mm}(num))

   Access matrix by number. For example::
     
     julia> mdlist(builtin(3))
     "chebspec"

.. function:: mdlist(builtin(num1:num2, ...))

   Access matrix by range and combinations. For example::

    julia> mdlist(builtin(1:4, 6, 10:15))
    11-element Array{String,1}:
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

We can define our own groups using the function ``setgroup!`` and
remove a defined group with ``deletegroup!``.

User definied groups may use arbitrary patterns to declare subsets of all available matrices,
and are not restricted to simple lists of alternative names.

See also::

    help?> MatrixDepot
    search: MatrixDepot matrixdepot

      julia MatrixDepot

      Give access to a wealth of sample and test matrices and accompanying data. A set of matrices is generated locally (with arguments controlling the special case). Another set is loaded from one of
      the publicly accessible matrix collections SuiteSparse Matrix Collection (formerly University of Florida Matrix Collection) and the Matrix Market Collection.

      Access is like

      using MatrixDepot
      
      A = matrixdepot("hilb", 10) # locally generated hilbert matrix dimensions (10,10)
      
      A = matrixdepot("HB/1138_bus")     # named matrix of the SuiteSparse Collection
      A = matrixdepot(sp(1))             # same matrix using numerical id
      A = matrixdepot("Harwell*/*/1138_bus") # matrix from the Matrix Market Collection 
      
      md = mdopen("*/bfly")   # named matrix with some extra data
      A = md.A
      co = md.coord
      tx = md("Gname_10.txt")
      
      md = mdopen("gravity", 10, false) # localy generated example with rhs and solution
      A = md.A
      b = md.b
      x = md.x

      commands:

      mdinfo, listdir, listgroups, matrixdepot, mdopen, listdata, mdlist,
      metasymbols, setgroup!, deletegroup!.

      selector patterns:

      strings, string-patterns (using "*", "?", "[]", "/", "**"), regular expressions: for names
      builtin(42), user(3,5), sp(10:11,6,2833), mm(1), mm(:): to access by integer id or all
      sp(pattern), mm(pattern) to access corresponding (alternative) matrix for other collection

      predicate patterns:

      isboolean, isinteger, isreal, iscomplex
      isgeneral, issymmetric, ishermitian, isskew
      isbuiltin, isuser, islocal, isremote, isloaded, isunloaded
      issvdok
      keyword(string expression), logical, hasdata(symbol), @pred(expression)
