.. Matrix Depot documentation master file, created by
   sphinx-quickstart on Tue Dec 23 00:46:02 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



Matrix Depot
------------

.. image:: logo1.png 

A test matrix collection for Julia.

.. toctree:: 
   :maxdepth: 2

   matrices.rst

Usage
-----

Every matrix in the collection is represented by a string
``matrix_name``, for example, the Cauchy matrix is represented by
``"cauchy"`` and the Hilbert matrix is represented by ``"hilb"``.

The properties of the matrices in the collection are also symbolized
by strings ``propertry_name``. For example, the class of the symmetric
matrices is symbolized by ``"symmetric"``.

* ``matrixdepot()`` returns a list of all the matrices in the
  collection.

* ``matrixdepot(matrix_name, p1, p2, ...)`` returns a matrix specified
  by the query string ``matrix_name``. ``p1, p2, ...`` are input
  parameters depending on ``matrix_name``.

* ``matrixdepot(matrix_name)`` returns the parameter options and the
  properties of ``matrix_name``.

* ``matrixdepot(prop1, prop2, ...)`` returns a list of matrices with
  the property ``prop1``, ``prop2``, etc.

We can define our own properties using the macro ``@addproperty`` and
remove a defined property using ``@rmproperty``.

* ``@addproperty property_name = ["matrix1", "matrix2", "matrix3"]``

* ``@rmproperty property_name``

   


