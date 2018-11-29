.. _interfaces:

Interface to Test Collections
=============================

The internal database is loaded automatically when using the module:

.. code::

    julia> using MatrixDepot
    include group.jl for user defined matrix generators
    verify download of index files...
    used remote site is https://sparse.tamu.edu/?per_page=All
    populating internal database...

Interface to the SuiteSparse Matrix Collection (formerly UFL collection)
------------------------------------------------------------------------

Use ``M = matrixdepot(NAME)`` or ``md = mdopen(NAME); M = md.A``, where ``NAME``
is ``collection_name + '/' + matrix_name``, to download a test matrix from the
`SuiteSparse Matrix Collection. https://sparse.tamu.edu/
For example::

  julia> md = mdopen("SNAP/web-Google")
  PG SNAP/web-Google(#2301)  916428x916428(5105039) 2002 [A] 'Directed Graph' [Web graph from Google]()

.. note:: 
   ``listnames("*/*")`` displays all the matrix names in the
   collection, including the newly downloaded matrices. All the matrix 
   data can be found by ``listnames("**")``.


If the matrix name is unique in the collections, we could also use
``matrixdepot(matrix_name)`` to download the data. If more than
one matrix has the same name, an error is thrown.

When download is complete, we can check matrix information using:

.. code::

  julia> mdinfo("SNAP/web-Google")
  SNAP/web-Google
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

  MatrixMarket matrix coordinate pattern general

  ──────────────────────────────────────────────────────────────────────

    •    UF Sparse Matrix Collection, Tim Davis

    •    http://www.cise.ufl.edu/research/sparse/matrices/SNAP/web-Google

    •    name: SNAP/web-Google

    •    [Web graph from Google]

    •    id: 2301

    •    date: 2002

    •    author: Google

    •    ed: J. Leskovec

    •    fields: name title A id date author ed kind notes

    •    kind: directed graph

  ───────────────────────────────────────────────────────────────────────
  ...


and generate it by accessing the field `A`.

.. code::

    julia> M = md.A
    916428×916428 SparseMatrixCSC{Bool,Int64} with 5105039 stored entries:
      [11343 ,      1]  =  true
      [11928 ,      1]  =  true
      [15902 ,      1]  =  true
      [29547 ,      1]  =  true
      [30282 ,      1]  =  true
      ⋮
      [788476, 916427]  =  true
      [822938, 916427]  =  true
      [833616, 916427]  =  true
      [417498, 916428]  =  true
      [843845, 916428]  =  true


You can convert the boolean pattern matrix to integer by ``M * 1``.

The metadata of a given matrix can be obtained by accessing properties of `md`.


Which properties are available is shown in the `md::MatrixDescriptor`:

.. code::

  julia> md = mdopen("TKK/t520")
  (IS TKK/t520(#1908)  5563x5563(286341/145952) 2008 [A, b, coord] 'Structural Problem' [T-beam, L = 520 mm, Quadratic four node DK type elements.  R Kouhia]()

and also by the special function ``metasymbols``:
  
.. code::

  julia> metasymbols(md)
  (:A, :b, :coord)

When you access a single matrix with ``matrixdepot(pattern)`` or ``mdopen(pattern)`` the full
matrix data are dowloaded implicitly in the background, if not yet available on the local disk
cache. 

When you access matrix information with ``mdinfo(pattern)`` for one or more matrices, the header
data of the matrix are downloaded implicitly, if not yet available on the local disk cache.

It is also possible to dowload a bulk of matrix data by ``MatrixDepot.loadinfo(pattern)`` and
``MatrixDepot.load(pattern)`` to populate the disk cache in advance of usage.


Interface to NIST Matrix Market
-------------------------------

Use ``M = matrixdepot(NAME)`` or ``md = mdopen(NAME); M = md.A``, where ``NAME``
is ``collection name + '/' + set name + '/' + matrix name`` to download a
test matrix from NIST Matrix Market:
http://math.nist.gov/MatrixMarket/. For example::

  julia> md = mdopen("Harwell-Boeing/lanpro/nos5")
  The collection-name and set-name may as always be replaced by wildcard patterns "*",
  as long as there exists only on name matching the pattern.

  julia> md = mdopen("*/*/bp__1400")
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
  100 28192  100 28192    0     0   4665      0  0:00:06  0:00:06 --:--:-- 10004
  download:/home/.../MatrixDepot/data/mm/Harwell-Boeing/smtape/bp__1400.mtx.gz

  (RG Harwell-Boeing/smtape/bp__1400(#M93)  822x822(4790)  [A] '' []()

 
Checking matrix information and generating matrix data are similar to 
the above case::

  julia> mdinfo(md) # or mdinfo("*/*/bp__1400")
    Harwell-Boeing/smtape/bp__1400
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    MatrixMarket matrix coordinate real general

    822 822 4790

There is no header information in this collection besides m, n, and dnz.

.. code::

  julia> md.A # or matrixdepot("Harwell-Boeing/smtape/bp__1400") 
  822x822 sparse matrix with 4790 Float64 entries:
	[1  ,   1]  =  1.0
	[1  ,   2]  =  0.001
	[26 ,   2]  =  -1.0
	[1  ,   3]  =  0.6885
	[25 ,   3]  =  0.9542
	[692,   3]  =  1.0
	[718,   3]  =  5.58
	⋮
	[202, 820]  =  -1.0
	[776, 820]  =  1.0
	[1  , 821]  =  0.4622
	[25 , 821]  =  0.725
	[28 , 821]  =  1.0
	[202, 821]  =  -1.0
	[796, 821]  =  1.0
	[2  , 822]  =  1.0

