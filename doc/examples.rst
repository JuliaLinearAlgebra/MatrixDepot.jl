.. _examples:

Examples
--------

To see all the matrices in the collection, type

.. code:: julia

    matrixdepot()

.. parsed-literal::

          | symmetric |  inverse  | ill-cond  |  pos-def  |   eigen   |
       kms|     *     |     *     |     *     |     *     |           |
      vand|           |     *     |     *     |           |           |
   poisson|     *     |     *     |           |     *     |     *     |
     frank|           |           |     *     |           |     *     |
     minij|     *     |     *     |           |     *     |     *     |
   clement|     *     |     *     |           |           |     *     |
   tridiag|     *     |     *     |     *     |     *     |     *     |
    circul|     *     |           |           |     *     |     *     |
  dingdong|     *     |           |           |           |     *     |
  hadamard|           |     *     |           |           |     *     |
     moler|     *     |     *     |     *     |     *     |           |
     invol|           |     *     |     *     |           |     *     |
   fiedler|     *     |     *     |           |           |     *     |
  binomial|           |           |           |           |           |
    lehmer|     *     |     *     |           |     *     |           |
    rohess|           |           |           |           |           |
   invhilb|     *     |     *     |     *     |     *     |           |
    lotkin|           |     *     |     *     |           |     *     |
 wilkinson|     *     |           |           |           |     *     |
      triw|           |     *     |     *     |           |           |
   randsvd|     *     |           |     *     |           |           |
    rosser|           |           |     *     |           |     *     |
     rando|           |           |           |           |           |
     magic|           |     *     |           |           |           |
     kahan|           |     *     |     *     |           |           |
    pascal|     *     |     *     |     *     |     *     |     *     |
  chebspec|           |           |           |           |     *     |
      chow|           |           |           |           |     *     |
  sampling|           |           |           |           |     *     |
      hilb|     *     |     *     |     *     |     *     |           |
    cauchy|     *     |     *     |     *     |     *     |           |
    wathen|     *     |           |           |     *     |     *     |
       pei|     *     |     *     |     *     |     *     |           |
    parter|           |           |           |           |     *     |
  forsythe|           |     *     |     *     |           |     *     |
  randcorr|     *     |           |           |           |           |
   neumann|           |           |           |           |     *     |
     grcar|           |           |           |           |     *     |


The meaning of the column heading is as follows:

-  ``"symmetric"``: the matrix is symmetric for some parameter values.

-  ``"inverse"``: the inverse of the matrix is known explicitly.

-  ``"ill-cond"``: the matrix is ill-conditioned for some parameter
   values.

-  ``"pos-def"``: the matrix is symmetric positive definite for some
   parameter values.

-  ``"eigen"``: the eigensystem of the matrix has some known results
   (explicit formulas for eigenvalues, eigenvectors, bounds of
   eigenvalues, etc).

We can generate a Hilbert matrix of size 4 by typing

.. code:: julia

    matrixdepot("hilb", 4)



.. parsed-literal::

    4x4 Array{Float64,2}:
     1.0       0.5       0.333333  0.25    
     0.5       0.333333  0.25      0.2     
     0.333333  0.25      0.2       0.166667
     0.25      0.2       0.166667  0.142857



and generate a circul matrix of size 5 by

.. code:: julia

    matrixdepot("circul", 5)



.. parsed-literal::

    5x5 Array{Float64,2}:
     1.0  2.0  3.0  4.0  5.0
     5.0  1.0  2.0  3.0  4.0
     4.0  5.0  1.0  2.0  3.0
     3.0  4.0  5.0  1.0  2.0
     2.0  3.0  4.0  5.0  1.0



We can type the matrix name to see the parameter options or matrix
properties.

.. code:: julia

    matrixdepot("hilb")

.. parsed-literal::

   Hilbert matrix: 
             
    Input options: 
             
    [type,] dim: the dimension of the matrix
             
    [type,] row_dim, col_dim: the row and column dimension 
             
    ['inverse', 'ill-cond', 'symmetric', 'pos-def']



.. code:: julia

    matrixdepot("hadamard")

.. parsed-literal::

   Hadamard matrix: 
             
    Input options: 
             
    [type,] dim: the dimension of the matrix, n is a power of 2 
             
    ['inverse', 'orthogonal', 'eigen']



From the information given, we notice that we can create a 4-by-6
rectanglular Hilbert matrix by

.. code:: julia

    matrixdepot("hilb", 4, 6)



.. parsed-literal::

    4x6 Array{Float64,2}:
     1.0       0.5       0.333333  0.25      0.2       0.166667
     0.5       0.333333  0.25      0.2       0.166667  0.142857
     0.333333  0.25      0.2       0.166667  0.142857  0.125   
     0.25      0.2       0.166667  0.142857  0.125     0.111111



We can aslo specify the data type

.. code:: julia

    matrixdepot("hilb", Float16, 5, 3)



.. parsed-literal::

    5x3 Array{Float16,2}:
     1.0      0.5      0.33325
     0.5      0.33325  0.25   
     0.33325  0.25     0.19995
     0.25     0.19995  0.16663
     0.19995  0.16663  0.14282



By inputing a matrix name, we can see what properties that matrix have.
Conversely, if we input a property (or properties), we can see all the
matrices (in the collection) having that property (or properties).

.. code:: julia

    matrixdepot("symmetric")



.. parsed-literal::

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

.. code:: julia

    matrixdepot("symmetric", "ill-cond")



.. parsed-literal::

    7-element Array{ASCIIString,1}:
     "hilb"   
     "cauchy" 
     "invhilb"
     "moler"  
     "pascal" 
     "pei"    
     "tridiag"



.. code:: julia

    matrixdepot("inverse", "ill-cond", "symmetric")



.. parsed-literal::

    7-element Array{ASCIIString,1}:
     "hilb"   
     "cauchy" 
     "invhilb"
     "moler"  
     "pascal" 
     "pei"    
     "tridiag"



Given a property, we can loop through all the matrices having this
propery

.. code:: julia

    # Multiply all matrices of the class "symmetric" and "ill-cond" and "inverse"
    A = eye(4)
    print("Identity matrix")
    for mat in intersect(matrixdepot("symmetric"), matrixdepot("ill-cond"), matrixdepot("inverse")) 
        print(" x $mat matrix")
        A = A * full(matrixdepot(mat, 4))    
    end
    println(" =")
    A    

.. parsed-literal::

    Identity matrix x hilb matrix x cauchy matrix x invhilb matrix x moler matrix x pascal matrix x pei matrix x tridiag matrix =




.. parsed-literal::

    4x4 Array{Float64,2}:
     153.12    -11.919    -15.4345   296.937
     109.896    -8.91857  -11.5976   214.433
      86.7524   -7.15714   -9.32857  169.702
      71.9139   -5.98707   -7.81497  140.876



The loop above can also be written as

.. code:: julia

    A = eye(4)
    print("Identity matrix")
    for mat in matrixdepot("symmetric", "ill-cond", "inverse")
        print(" x $mat matrix")
        A = A * full(matrixdepot(mat, 4))
    end
    println(" =")
    A

.. parsed-literal::

    Identity matrix x hilb matrix x cauchy matrix x invhilb matrix x moler matrix x pascal matrix x pei matrix x tridiag matrix =




.. parsed-literal::

    4x4 Array{Float64,2}:
     153.12    -11.919    -15.4345   296.937
     109.896    -8.91857  -11.5976   214.433
      86.7524   -7.15714   -9.32857  169.702
      71.9139   -5.98707   -7.81497  140.876



User Defined Properties
-----------------------

We can define properties in MatrixDepot. Since each property in Matrix
Depot is a list of strings, you can simply do, for example,

.. code:: julia

    spd = matrixdepot("symmetric", "pos-def")



.. parsed-literal::

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



.. code:: julia

    myprop = ["lehmer", "cauchy", "hilb"]



.. parsed-literal::

    3-element Array{ASCIIString,1}:
     "lehmer"
     "cauchy"
     "hilb"  



Then use it in your tests like

.. code:: julia

    for matrix in myprop
        A = matrixdepot(matrix, 6)
        L, U, p = lu(A) #LU factorization
        err = norm(A[p,:] - L*U, 1) # 1-norm error
        println("1-norm error for $matrix matrix is ", err)
    end    

.. parsed-literal::

    1-norm error for lehmer matrix is 1.1102230246251565e-16
    1-norm error for cauchy matrix is 5.551115123125783e-17
    1-norm error for hilb matrix is 2.7755575615628914e-17


To add a property permanently for future use, we put the macro
``@addproperty`` at the beginning.

.. code:: julia

    @addproperty myfav = ["lehmer", "cauchy", "hilb"]



.. parsed-literal::

    87



.. code:: julia

    @addproperty spd = matrixdepot("symmetric", "pos-def")



.. parsed-literal::

    195



We need to **restart** Julia to see the changes. Type

.. code:: julia

    matrixdepot()

.. parsed-literal::

    
              | symmetric |  inverse  | ill-cond  |  pos-def  |  eigen    |
          vand|           |     *     |     *     |           |           |
         frank|           |           |     *     |           |     *     |
         minij|     *     |     *     |           |     *     |     *     |
       clement|     *     |     *     |           |           |     *     |
       tridiag|     *     |     *     |     *     |     *     |     *     |
        circul|     *     |           |           |     *     |     *     |
      dingdong|     *     |           |           |           |     *     |
      hadamard|           |     *     |           |           |     *     |
         moler|     *     |     *     |     *     |     *     |           |
         invol|           |     *     |     *     |           |     *     |
       fiedler|     *     |     *     |           |           |     *     |
      binomial|           |           |           |           |           |
        lehmer|     *     |     *     |           |     *     |           |
       invhilb|     *     |     *     |     *     |     *     |           |
        lotkin|           |     *     |     *     |           |     *     |
          triw|           |     *     |     *     |           |           |
         magic|           |     *     |           |           |           |
         kahan|           |     *     |     *     |           |           |
        pascal|     *     |     *     |     *     |     *     |     *     |
      chebspec|           |           |           |           |     *     |
          hilb|     *     |     *     |     *     |     *     |           |
        cauchy|     *     |     *     |     *     |     *     |           |
           pei|     *     |     *     |     *     |     *     |           |
      forsythe|           |     *     |     *     |           |     *     |
         grcar|           |           |           |           |     *     |
    
    New Properties:
    
    spd = [ hilb, cauchy, circul, invhilb, moler, pascal, pei, minij, tridiag, lehmer, ] 
    
    myfav = [ lehmer, cauchy, hilb, ] 
    


Notice new defined properties have been included. We can use them as

.. code:: julia

    matrixdepot("myfav")



.. parsed-literal::

    3-element Array{ASCIIString,1}:
     "lehmer"
     "cauchy"
     "hilb"  



We can remove a property using the macro ``@rmproperty``. As before, we
need to **restart** Julia to see the changes.

.. code:: julia

    @rmproperty myfav



.. parsed-literal::

    153



.. code:: julia

    matrixdepot()

.. parsed-literal::

    
              | symmetric |  inverse  | ill-cond  |  pos-def  |  eigen    |
          vand|           |     *     |     *     |           |           |
         frank|           |           |     *     |           |     *     |
         minij|     *     |     *     |           |     *     |     *     |
       clement|     *     |     *     |           |           |     *     |
       tridiag|     *     |     *     |     *     |     *     |     *     |
        circul|     *     |           |           |     *     |     *     |
      dingdong|     *     |           |           |           |     *     |
      hadamard|           |     *     |           |           |     *     |
         moler|     *     |     *     |     *     |     *     |           |
         invol|           |     *     |     *     |           |     *     |
       fiedler|     *     |     *     |           |           |     *     |
      binomial|           |           |           |           |           |
        lehmer|     *     |     *     |           |     *     |           |
       invhilb|     *     |     *     |     *     |     *     |           |
        lotkin|           |     *     |     *     |           |     *     |
          triw|           |     *     |     *     |           |           |
         magic|           |     *     |           |           |           |
         kahan|           |     *     |     *     |           |           |
        pascal|     *     |     *     |     *     |     *     |     *     |
      chebspec|           |           |           |           |     *     |
          hilb|     *     |     *     |     *     |     *     |           |
        cauchy|     *     |     *     |     *     |     *     |           |
           pei|     *     |     *     |     *     |     *     |           |
      forsythe|           |     *     |     *     |           |     *     |
         grcar|           |           |           |           |     *     |
    
    New Properties:
    
    spd = [ hilb, cauchy, circul, invhilb, moler, pascal, pei, minij, tridiag, lehmer, ] 
    


More Examples
-------------

An interesting test matrix is magic square. It can be generated as

.. code:: julia

    M = matrixdepot("magic", 5)



.. parsed-literal::

    5x5 Array{Int64,2}:
     17  24   1   8  15
     23   5   7  14  16
      4   6  13  20  22
     10  12  19  21   3
     11  18  25   2   9



.. code:: julia

    sum(M,1)



.. parsed-literal::

    1x5 Array{Int64,2}:
     65  65  65  65  65



.. code:: julia

    sum(M,2)



.. parsed-literal::

    5x1 Array{Int64,2}:
     65
     65
     65
     65
     65



.. code:: julia

    sum(diag(M))



.. parsed-literal::

    65



.. code:: julia

    p = [5:-1:1]
    sum(diag(M[:,p]))



.. parsed-literal::

    65



Pascal Matrix can be generated as

.. code:: julia

    P = matrixdepot("pascal", 6)



.. parsed-literal::

    6x6 Array{Int64,2}:
     1  1   1   1    1    1
     1  2   3   4    5    6
     1  3   6  10   15   21
     1  4  10  20   35   56
     1  5  15  35   70  126
     1  6  21  56  126  252



Notice the Cholesky factor of the Pascal matrix has Pascal's triangle
rows.

.. code:: julia

    chol(P)



.. parsed-literal::

    6x6 Array{Float64,2}:
     1.0  1.0  1.0  1.0  1.0   1.0
     0.0  1.0  2.0  3.0  4.0   5.0
     0.0  0.0  1.0  3.0  6.0  10.0
     0.0  0.0  0.0  1.0  4.0  10.0
     0.0  0.0  0.0  0.0  1.0   5.0
     0.0  0.0  0.0  0.0  0.0   1.0


