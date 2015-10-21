.. _user:

Adding New Matrix Generators
============================

Matrix Depot provides a diverse collection of 
test matrices, including parametrized matrices
and real-life matrices. But occasionally, you 
may want to define your own matrix generators and 
be able to use them from Matrix Depot. Can I do that? 
Yes and it is very simple :)

Declaring Generators
--------------------

To get a feel of how it works, let's see an example. 
Suppose we want to include a random symmetric matrix ``randsym``
in Matrix Depot::

  function randsym{T}(::Type{T}, n)
    A = zeros(T, n, n)
    for j = 1:n
        for i = j:n
            A[i,j] = randn()
        end
    end
    A = A + tril(A, -1)'
    return A
  end

We first need to find out where Matrix Depot is installed. This 
can be done by::

  julia> Pkg.dir("MatrixDepot")
  "/home/weijian/.julia/v0.4/MatrixDepot"

For me, the package is installed at
``/home/weijian/.julia/v0.4/MatrixDepot``. Now we open the file
``user/user.jl`` inside ``MatrixDepot``. It looks like this::

  usermatrixclass =
  @compat Dict(



  );

  ##########################################
  # Please put your matrix generators here #
  ##########################################

We can copy and paste the function ``randsym`` anywhere below the
comments and use the function ``include_generator`` to include it::
  
  usermatrixclass =
  @compat Dict(



  );

  ##########################################
  # Please put your matrix generators here #
  ##########################################

  function randsym{T}(::Type{T}, n)
   A = zeros(T, n, n)
    for j = 1:n
        for i = j:n
            A[i,j] = randn()
        end
    end
    A = A + tril(A, -1)'
    return A
  end
  include_generator(FunctionName, "randsym", randsym)

This is it. We can now use it from Matrix Depot::

  julia> matrixdepot()

  Matrices:
    1) baart            2) binomial         3) blur             4) cauchy        
    5) chebspec         6) chow             7) circul           8) clement       
    9) deriv2          10) dingdong        11) fiedler         12) forsythe      
   13) foxgood         14) frank           15) golub           16) gravity       
   17) grcar           18) hadamard        19) hankel          20) heat          
   21) hilb            22) invhilb         23) invol           24) kahan         
   25) kms             26) lehmer          27) lotkin          28) magic         
   29) minij           30) moler           31) neumann         32) oscillate     
   33) parter          34) pascal          35) pei             36) phillips      
   37) poisson         38) prolate         39) randcorr        40) rando         
   41) randsvd         42) randsym         43) rohess          44) rosser        
   45) sampling        46) shaw            47) spikes          48) toeplitz      
   49) tridiag         50) triw            51) vand            52) wathen        
   53) wilkinson       54) wing
  Groups:
    all           data          eigen         ill-cond    
    inverse       pos-def       random        regprob     
    sparse        symmetric 

  julia> matrixdepot("randsym", 5)
  5x5 Array{Float64,2}:
   1.57579    0.474591  0.0261732  -0.536217  -0.0900839
   0.474591   0.388406  0.77178     0.239696   0.302637 
   0.0261732  0.77178   1.7336      1.72549    0.127008 
  -0.536217   0.239696  1.72549     0.304016   1.5854   
  -0.0900839  0.302637  0.127008    1.5854    -0.656608 

  julia> matrixdepot("randsym", Float32, 5)
  5x5 Array{Float32,2}:
  -0.633797  -0.154157   0.972601  0.554571  -0.692858
  -0.154157  -0.319152  -0.710942  2.81623    1.2637  
   0.972601  -0.710942  -0.165526  1.16547   -0.705227
   0.554571   2.81623    1.16547   0.351268   0.410586
  -0.692858   1.2637    -0.705227  0.410586  -0.786438

To make it more useful, we can include the helper strings and group information::

  helplines = "random symmetric matrix:
            \n Input options: [type, n]: the dimension of the matrix is n."
  include_generator(Help, helplines, randsym)
  include_generator(Group, "random", randsym)
  include_generator(Group, "symmetric", randsym)

Now we can do::

  julia> matrixdepot("randsym")
  random symmetric matrix:
            
  Input options: [type, n]: the dimension of the matrix is n.

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



