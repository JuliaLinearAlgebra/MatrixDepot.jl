Usage
-----

Every matrix in the collection is represented by a string
``matrix_name``, for example, the Cauchy matrix is represented by
``"cauchy"`` and the Hilbert matrix is represented by ``"hilb"``.

The properties of the matrices in the collection are also symbolized
by strings ``propertry_name``. For example, the class of the symmetric
matrices is symbolized by ``"symmetric"``.

.. function:: matrixdepot()

  Return a list of all the matrices in the collection.

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

  Return the parameter options and the properties of ``matrix_name``. 
  For example::

    julia> matrixdepot("moler")
    Moler Matrix:
             
    Input options:
             
    (type), dim, alpha: dim is the dimension of the matrix,
                alpha is a scalar.
             
    (type), dim: alpha = -1.
             
    ['inverse', 'ill-cond', 'symmetric', 'pos-def']


.. function:: matrixdepot(group_name)

  Return a list of matrices with the property `group_name`. For
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

.. function:: matrixdepot(prop1, prop2, ...)

  Return a list of matrices with the property ``prop1``, ``prop2``, etc. 
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

   Access matrix by ``UnitRange{Int64}``. For example::

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

    Read the matrix data given by ``name``.

We can define our own groups using the macro ``@addgroup`` and
remove a defined group using ``@rmgroup``.

.. function:: @addgroup group_name = ["matrix1", "matrix2", "matrix3"]

  Create a new group `"group_name"` such that `"matrix1"`, `"matrix2"`
  and `"matrix3"` have this group.

.. function:: @rmgroup group_name
  
   Delete a created group `group_name`.
