.. Matrix Depot documentation master file, created by
   sphinx-quickstart on Tue Dec 23 00:46:02 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. toctree:: 
   :hidden:

   matrices.rst
   examples.rst
   properties.rst
   interface.rst

MatrixDepot.jl
==============

**Contents:**

* :ref:`install` 
* :ref:`usage` 
* :ref:`matrices` 
* :ref:`properties`
* :ref:`examples` 
* :ref:`interfaces`

MatrixDepot.jl is an open source test matrix collection for Julia.
Source and license information can be found at `GitHub`_.

.. _GitHub: https://github.com/weijianzhang/MatrixDepot.jl

.. _install:

Install
-------

To install the release version, type::

  julia> Pkg.add("MatrixDepot")

To install the latest development version, type::

  julia> Pkg.clone("MatrixDepot")

.. _usage:

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


.. function:: matrixdepot(property_name)

  Return a list of matrices with the property `property_name`. For
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

.. function:: MatrixDepot.get(name [,newname] [, collection =:UF])

    Download a matrix from test matrix collections, where
    ``name`` is a string of collection name + ``/`` + matrix name if
    ``collection = :UF`` (UF Sparse Matrix Collection) and ``name`` is
    collection name + ``/`` + set name + ``/`` + matrix name if
    ``collection = :MM`` (Matrix Market). 
    If ``newname`` is present, the downloaded matrix will be renamed to 
    ``newname``.
    For example::
      
      julia> MatrixDepot.get("HB/illc1850")
      julia> MatrixDepot.get("Harwell-Boeing/lanpro/nos5", collection = :MM)
      

.. function:: MatrixDepot.update()

    Update matrix collection database from the web server.


.. function:: matrixdepot(name)

    Output matrix information, where ``name`` is a matrix data.

.. function:: matrixdepot(name, :r)

    Read the matrix data given by ``name``.

We can define our own properties using the macro ``@addproperty`` and
remove a defined property using ``@rmproperty``.

.. function:: @addpropery property_name = ["matrix1", "matrix2", "matrix3"]

  Create a new property `"property_name"` such that `"matrix1"`, `"matrix2"`
  and `"matrix3"` have this property.

.. function:: @rmproperty property_name
  
   Delete a created property `property_name`.



   


