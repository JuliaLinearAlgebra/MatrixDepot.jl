
.. _properties:

Groups
======

Groups are lists of matrix names and we use them to
categorize matrices in Matrix Depot. The list below shows
all the predefined groups in Matrix Depot and we can extend
this list by defining new groups. Group names are noted as 
symbols, e.g. `:symmetric`.

Predefined Groups
-----------------

.. glossary::
   :sorted:

   symmetric
      The matrix is symmetric for some parameter values.

   inverse
      The inverse of the matrix is known explicitly.

   illcond
       The matrix is ill-conditioned for some parameter values.

   posdef
       The matrix is positive definite for some parameter values.

   eigen
       Part of the eigensystem of the matrix is explicitly known.

   sparse
      The matrix is sparse.

   random
      The matrix has random entries.

   data
      The matrix has been downloaded from UF sparse matrix collection or
      the Matrix Market collection. 

   regprob
      The output is a test problem for Regularization Methods.

   all
      All the matrices in the collection. 

   graph
      An adjacency matrix of a graph.  

Adding New Groups
-----------------

New groups can be added with the macro ``@addgroup``in the REPL::

    @addgroup myfav = ["lehmer", "cauchy", "hilb"]
    @addgroup test_for_paper2 = ["tridiag", "sampling", "wing"]

    julia> listgroups()
    14-element Array{Symbol,1}:
    :all
    :builtin
    :local
    :user
    :eigen
    :graph
    :illcond
    :inverse
    :posdef
    :random
    :regprob
    :sparse
    :symmetric
    :myfav

    julia> listnames(:myfav)
    list(3)            
    ––––––– –––– ––––––
    cauchy  hilb lehmer

