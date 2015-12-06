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

When Matrix Depot is first loaded, a new directory ``myMatrixDepot``
will be created. Matrix Depot automatically includes all Julia files
in this directory. Hence, all we need to do is to copy
the generator files to ``path/to/MatrixDepot/myMatrixDepot`` and use
the function ``include_generator`` to declare them.

.. function:: include_generator(Stuff To Be Included, Stuff, f)

   Includes a piece of information of the function ``f`` to Matrix Depot,
   where ``Stuff To Be Included`` is one of the following:
   
    * ``FunctionName``: the function name of ``f``. In this case, 
      ``Stuff`` is a string representing ``f``.
 
    * ``Group``: the group where ``f`` belongs. In this case, 
      ``Stuff`` is the group name.

Examples
--------- 

To get a feel of how it works, let's see an example. 
Suppose we have a file ``myrand.jl`` which contains two 
matrix generator ``randsym`` and ``randorth``::

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
        end
      end
      A = A + tril(A, -1)'
      return A
   end

   """
   random Orthogonal matrix
   ========================

   *Input options:*

   + n: the dimension of the matrix
     """	
     randorth(n) = qr(randn(n,n))[1]

We first need to find out where Matrix Depot is installed. This 
can be done by::

  julia> Pkg.dir("MatrixDepot")
  "/home/weijian/.julia/v0.4/MatrixDepot"

For me, the package is installed at
``/home/weijian/.julia/v0.4/MatrixDepot``. We can copy ``myrand.jl``
to ``/home/weijian/.julia/v0.4/MatrixDepot/myMatrixDepot``. 
Now we open the file
``myMatrixDepot/generator.jl`` and write::

  include_generator(FunctionName, "randsym", randsym)
  include_generator(FunctionName, "randorth", randorth)

This is it. We can now use them from Matrix Depot::

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
  41) rando           42) randorth        43) randsvd         44) randsym       
  45) rohess          46) rosser          47) sampling        48) shaw          
  49) spikes          50) toeplitz        51) tridiag         52) triw          
  53) ursell          54) vand            55) wathen          56) wilkinson     
  57) wing          
 Groups:
  all           data          eigen         ill-cond    
  inverse       pos-def       random        regprob     
  sparse        symmetric  

  julia> matrixdepot("randsym")
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

We can also add group information::

  include_generator(Group, "random", randsym)
  include_generator(Group, "symmetric", randsym)

Now if we type::

  julia> matrixdepot("random")
  9-element Array{ASCIIString,1}:
  "golub"    
  "oscillate"
  "randcorr" 
  "rando"    
  "randsvd"  
  "randsym"  
  "rohess"   
  "rosser"   
  "wathen" 

  julia> matrixdepot("symmetric")
  22-element Array{ASCIIString,1}:
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
  ⋮          
  "pascal"   
  "pei"      
  "poisson"  
  "prolate"  
  "randcorr" 
  "randsym"  
  "tridiag"  
  "wathen"   
  "wilkinson"

the function ``randsym`` is now part of the group ``symmetric`` and
``random``.


It is a good idea to back up your changes. For example, we 
could save it on GitHub by creating a new repository named ``myMatrixDepot``.
(See https://help.github.com/articles/create-a-repo/ for details of creating a new repository on GitHub.)
Then we go to the directory ``path/to/MatrixDepot/myMatrixDepot`` and type::

  git init
  git add group.jl
  git add generator.jl
  git commit -m "first commit"
  git remote add origin https://github.com/your-user-name/myMatrixDepot.git
  git push -u origin master

  


