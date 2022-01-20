.. _interfaces:

Interface to Test Collections
=============================

The internal database is loaded automatically when using the module::

  julia> using MatrixDepot
  [ Info: verify download of index files...
  [ Info: reading database
  [ Info: adding metadata...
  [ Info: adding svd data...
  [ Info: writing database
  [ Info: used remote sites are sparse.tamu.edu with MAT index and math.nist.gov with HTML index

Interface to the SuiteSparse Matrix Collection (formerly UFL collection)
------------------------------------------------------------------------

Use ``M = matrixdepot(NAME)`` or ``md = mdopen(NAME); M = md.A``, where ``NAME``
is ``collection_name + '/' + matrix_name``, to download a test matrix from the
SuiteSparse Matrix Collection:
https://sparse.tamu.edu/

For example::

   julia> md = mdopen("SNAP/web-Google")
   PG SNAP/web-Google(#2301)  916428x916428(5105039) 2002 [A] 'Directed Graph' [Web graph from Google]()

.. note:: 
   ``listnames("*/*")`` displays all the matrix names in the
   collection, including the newly downloaded matrices. All the matrix 
   data can be found by ``listnames("**")``.


If the matrix name is unique in the collections, we could also use
``matrixdepot("**/matrix_name")`` to download the data. If more than
one matrix has the same name, an error is thrown.

When download is complete, we can check matrix information using:

.. code::

  julia> mdinfo(md) # - or mdinfo("SNAP/web-Google")
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


and generate it by accessing the field ``A``.

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

The metadata of a given matrix can be obtained by accessing properties of ``md``
or ``md.data``.


Which properties are available is shown in the ``MatrixDescriptor``:

.. code::

  julia> md = mdopen("TKK/t520")
  (IS TKK/t520(#1908)  5563x5563(286341/145952) 2008 [A, b, coord] 'Structural Problem' [T-beam, L = 520 mm, Quadratic four node DK type elements.  R Kouhia]()

and also by the special function ``metasymbols``:
  
.. code::

    julia> metasymbols(md)
    (:A, :b, :coord)

You can obtain the metatdata ``A``, ``b``, and ``coord`` by::

    julia> md.coord
    5563×3 Array{Float64,2}:
    0.0   0.0   0.0
    0.0   0.0   0.0
    ⋮
    520.0  50.0  15.0
    520.0  50.0  15.0


When you access a single matrix with ``matrixdepot(pattern)`` or ``mdopen(pattern)``
the full matrix data are dowloaded implicitly in the background, if not yet available
on the local disk cache. 

When you access matrix information with ``mdinfo(pattern)`` for one or more matrices,
the header data of the matrix are downloaded implicitly, if not yet available on the
local disk cache.

It is also possible to dowload a bulk of matrix data by ``MatrixDepot.loadinfo(pattern)``
and ``MatrixDepot.load(pattern)`` to populate the disk cache in advance of usage.
If you want to access the Singular Value Decomposition (svd) data available for quite a
few of the Suite Sparse collection, you explicitly have to use
``MatrixDepot.loadsvd(pattern)``.

The following example demonstrates how to access SVD data (derived from singular
value decomposition of the matrix). The predicate ``issvdok`` selects all
matrices which have SVD data loaded. 

.. code::

    julia> mdlist(issvdok & @pred(5900 <= n < 6000))
    8-element Array{String,1}:
    "AG-Monien/ukerbe1"
    "Cote/mplate"
    "HB/man_5976"
    "Hamrle/Hamrle2"
    "JGD_Homology/cis-n4c6-b4"
    "JGD_Homology/n4c6-b4"
    "Schenk_IBMNA/c-32"
    "TOKAMAK/utm5940"

    julia> md = mdopen("AG-Monien/ukerbe1")
    (PS AG-Monien/ukerbe1(#2422)  5981x5981(15704/7852) 1998 [A, coord] '2D/3D problem' [2D finite element problem])()

    julia> reshape(propertynames(md.data),11, 4)
    11×4 Array{Symbol,2}:
    :name      :ed                  :amd_vnz             :xmax
    :id        :fields              :amd_rnz             :svdok
    :metadata  :notes               :amd_flops           :norm
    :m         :nnzdiag             :ncc                 :minsv
    :n         :pattern_symmetry    :nblocks             :cond
    :nnz       :numerical_symmetry  :sprank              :rank
    :dnz       :posdef              :lowerbandwidth      :nullspace
    :kind      :isND                :upperbandwidth      :svgap
    :date      :isGraph             :rcm_lowerbandwidth  :svdstatus
    :title     :cholcand            :rcm_upperbandwidth  :svdhow
    :author    :amd_lnz             :xmin                :sv

    julia> md.data.rank
    4108

    julia> sv = md.data.sv
    5981-element Array{Float64,1}:
    3.131396665809681
    3.1313966657795302
    3.0773931783051283
    ⋮
    3.936758260137112e-18
    1.550044427797539e-18
    7.503077983559783e-19
    8.317116401465794e-22

For the meaning of the property names see also: https://sparse.tamu.edu/statistics.

Interface to NIST Matrix Market
-------------------------------

Use ``M = matrixdepot(NAME)`` or ``md = mdopen(NAME); M = md.A``, where ``NAME``
is ``collection name + '/' + set name + '/' + matrix name`` to download a
test matrix from NIST Matrix Market:
http://math.nist.gov/MatrixMarket/. For example::

  julia> md = mdopen("Harwell-Boeing/lanpro/nos5")

The collection-name and set-name may as always be replaced by wildcard patterns "*",
as long as there exists only one name matching the pattern.

.. code::

  julia> md = mdopen("*/*/bp__1400")
  download:/home/.../MatrixDepot/data/mm/Harwell-Boeing/smtape/bp__1400.mtx.gz

  (RG Harwell-Boeing/smtape/bp__1400(#M93)  822x822(4790)  [A] '' []()

 
Checking matrix information and generating matrix data are similar to 
the above case::

  julia> mdinfo(md) # or mdinfo("*/*/bp__1400")
    Harwell-Boeing/smtape/bp__1400
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    MatrixMarket matrix coordinate real general

    822 822 4790

There is typically no header information in this collection besides m, n, and dnz.

.. code::

  julia> md.A # or matrixdepot("Harwell-Boeing/smtape/bp__1400") 
  822x822 sparse matrix with 4790 Float64 entries:
	[1  ,   1]  =  1.0
	[1  ,   2]  =  0.001
	[26 ,   2]  =  -1.0
	[1  ,   3]  =  0.6885
	⋮
	[25 , 821]  =  0.725
	[28 , 821]  =  1.0
	[202, 821]  =  -1.0
	[796, 821]  =  1.0
	[2  , 822]  =  1.0


Matrix Identification (Patterns)
--------------------------------

A ``pattern`` is used to select from the available problems. There are several elementary
and combination forms of patterns. Each pattern selects an array of matrix names
currently found in the database. Matrix names contain zero, one or two ``/`` characters.

Patterns can be used to select subsets or individual problems from the depot.
They are used as arguments of the functions::

    matrixdepot
    mdinfo
    mdlist
    listnames
    listdata
    listdir     - only strings
    mdopen      - single match required    
    MatrixDepot.loadinfo
    MatrixDepot.loadsvd

1. ``AbstractString`` with wildcard characters ``*``, ``**``, and ``?``.

  Here ``?`` stands for one arbitrary character in the matrix name excluding ``/``, and ``*``
  for a sequence of arbitray characters excluding ``/``. ``**`` stands for an arbitrary
  sequence of characters including ``/``.

  Example: ``"*/???/w*"``, ``"**/1138*"``

2. One of the integer identifiers ``builtin(n)``, ``user(n)``, ``sp(n)``, and ``mm(n)``.
  
  Here the respectively built-in, user-defined, suite-sparse, matrix-market problems are
  numbered. ``n`` may be a positive integer, a range of integers, or a list of the previous.

  Example: ``builtin(1,3,10:11)``

2. One of the patterns for accessing alternate collection ``mm(pattern), sp(pattern)``.

   If ``pattern`` selects a matrix form the Suite Sparse, then ``mm(pattern)` selects the
   corresponding matrix form Matrix Market. If no such matrix exists, nothing is returned.
   If the matrix selected by ``pattern`` is already in Matrix Market, the same matrix is returned.

   Example: ``sp("*/*/1138_bus") == ["HB/1138_bus"]``.

   Note:
   The matrix names may be mangled sometimes.

2. A ``Symbol`` indicating one of the defined groups.

  Example: ``:symmetric``

3. A ``Function`` (predicate of subtypes of ``MatrixData``).

  Example: ``data::MatrixDepot.RemoteMatrixData -> data.n <= 100``


4. A ``@pred`` predicate function.

  Example: ``@pred(n <= 100)`` which is a shorthand for the previous example.

5. One of the predefined predicate functions.

  ``isreal, iscomplex, isinteger, isboolean,``
  ``islocal, isremote,``
  ``isbuiltin, isuser, isloaded, isunloaded,``
  ``isgeneral, issymmetric, ishermitian, isskew,``
  ``issvdok, isposdef``

6. A list ``AbstractVector`` or disjunction of any number of all forms, meaning ``OR``.

  Example: ``[mm(1), sp(1)]`` or equivalently: ``mm(1) | sp(1)`` (note: single ``|``)

7. A tuple or conjuction of any number of all forms, meaning ``AND``.

  Example: ``(mm(:), @pred(m == 1000))`` or shorter: ``m(::) & @pred(m == 1000)``

8. The negation of any of the previous by a unary ``~`` and parenthesized terms

  Example: ``(issymmetric | ishermitian) & ~isposdef``

