.. _user:

Adding New Matrix Generators
============================

Matrix Depot provides a diverse collection of 
test matrices, including parametrized matrices
and real-life matrices. But occasionally, you 
may want to define your own matrix generators and 
be able to use them from Matrix Depot. 

Declaring Generators
--------------------

Users may override the ``include_generator`` function to declare new
generators. Note that ``include_generator`` must be called at runtime, so
users creating packages which define custom matrix test suites should place
their calls to ``include_generator`` in the ``__init__`` function of their
module to avoid
`precompilation issues <https://docs.julialang.org/en/v2/manual/modules/#Module-initialization-and-precompilation>`_.

.. function:: include_generator(Stuff_To_Be_Included, Stuff, f)

   Includes a piece of information of the function ``f`` to Matrix Depot,
   where ``Stuff_To_Be_Included`` is one of the following:
   
    * ``FunctionName``: the function name of ``f``. In this case, 
      ``Stuff`` is a string representing ``f``.
 
    * ``Group``: the group where ``f`` belongs. In this case, 
      ``Stuff`` is the group name.

Examples
--------- 

To get a feel of how it works, let's see an example. Suppose we have two
matrix generators ``randsym`` and ``randorth``::

  """
  random symmetric matrix
  =======================

  *Input options:* 

  + n: the dimension of the matrix
  """
  function randsym(n)
      A = zeros(n, n)
      for j = 1:n
        for i = j:n
          A[i,j] = randn()
          if i != j; A[j,i] = A[i,j] end
        end
      end
      return A
  end

  """
  random Orthogonal matrix
  ========================

  *Input options:*

  + n: the dimension of the matrix
  """
  randorth(n) = Matrix(qr(randn(n,n)).Q)

We can then write::

  MatrixDepot.include_generator(FunctionName, "randsym", randsym)
  MatrixDepot.include_generator(FunctionName, "randorth", randorth)

and when we are done including generators, we need to reinitializie the database::

  MatrixDepot.init()

This is it. We can now use them from Matrix Depot::

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
    –––––––––– –––––––––
    1 randorth 2 randsym

    Groups                                                          
    ––––––– ––––– ––––– ––––––– –––––– ––––––– –––––––––––––––      
    all     local eigen illcond posdef regprob symmetric            
    builtin user  graph inverse random sparse  test_for_paper2      

    Suite Sparse of  
    –––––––––––– ––––
    2773         2833

    MatrixMarket of 
    –––––––––––– –––
    488          498


    julia> mdinfo("randsym")
        random symmetric matrix
        ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

        Input options: 

        •  n: the dimension of the matrix

    julia> matrixdepot("randsym", 5)
    5x5 Array{Float64,2}:
    1.57579    0.474591  0.0261732  -0.536217  -0.0900839
    0.474591   0.388406  0.77178     0.239696   0.302637 
    0.0261732  0.77178   1.7336      1.72549    0.127008 
    -0.536217   0.239696  1.72549     0.304016   1.5854   
    -0.0900839  0.302637  0.127008    1.5854    -0.656608 

    julia> A = matrixdepot("randorth", 5)
    5x5 Array{Float64,2}:
    -0.359134   0.401435   0.491005  -0.310518   0.610218
    -0.524132  -0.474053  -0.53949   -0.390514   0.238764
    0.627656   0.223519  -0.483424  -0.104706   0.558054
    -0.171077   0.686038  -0.356957  -0.394757  -0.465654
    0.416039  -0.305802   0.326723  -0.764383  -0.205834

    julia> A'*A
    5x5 Array{Float64,2}:
    1.0           8.32667e-17   1.11022e-16   5.55112e-17  -6.93889e-17
    8.32667e-17   1.0          -1.80411e-16  -2.77556e-17  -5.55112e-17
    1.11022e-16  -1.80411e-16   1.0           1.94289e-16  -1.66533e-16
    5.55112e-17  -2.77556e-17   1.94289e-16   1.0           1.38778e-16
    -6.93889e-17  -5.55112e-17  -1.66533e-16   1.38778e-16   1.0 

We can also add group information with::

    MatrixDepot.include_generator(MatrixDepot.Group, :random, randsym)
    MatrixDepot.include_generator(MatrixDepot.Group, :symmetric, randsym)
    MatrixDepot.include_generator(MatrixDepot.Group, :random, randorth)
    MatrixDepot.init()

For example::

    julia> mdlist(:symmetric)
    22-element Array{String,1}:
    "cauchy"
    "circul"
    "clement"
    "dingdong"
    "fiedler"
    "hankel"
    "hilb"
    "invhilb"
    "kms"
    "lehmer"
    "minij"
    "moler"
    "oscillate"
    "pascal"
    "pei"
    "poisson"
    "prolate"
    "randcorr"
    "randsym"
    "tridiag"
    "wathen"
    "wilkinson"

    julia> listnames(:random)
    list(13)                                                           
    –––––––– ––––––––– –––––––– –––––––– ––––––– –––––––––– ––––––     
    erdrey   golub     randcorr randorth randsym rosser     wathen     
    gilbert  oscillate rando    randsvd  rohess  smallworld            

the function ``randsym`` will be part of the groups ``:symmetric`` and ``:random``
while ``randorth`` is in group ``:random``.