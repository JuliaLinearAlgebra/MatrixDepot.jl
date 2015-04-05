.. _interfaces:

Interface to Test Matrix Collections
====================================

Interface to the UF Sparse Matrix Collection
---------------------------------------------

Before downloading test matrices, we should first update the database::

  julia> MatrixDepot.update()

Use ``MatrixDepot.get(NAME)``, where ``NAME`` is ``collection name
+'/' + matrix name``,  to download a test matrix from the
`UF Sparse Matrix Collection <http://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html>`_.
For example::

  julia> MatrixDepot.get("SNAP/web-Google")

.. note:: 
   ``matrixdepot()`` displays all the matrices in the
   collection, including the newly downloaded matrices. All the matrix 
   data can be found by ``matrixdepot("data")``. 
	  

When download is complete, we can check matrix information using::

  julia> matrixdepot("SNAP/web-Google")

  %%MatrixMarket matrix coordinate pattern general
  %-------------------------------------------------------------------------------
  % UF Sparse Matrix Collection, Tim Davis
  % http://www.cise.ufl.edu/research/sparse/matrices/SNAP/web-Google
  % name: SNAP/web-Google
  % [Web graph from Google]
  % id: 2301
  % date: 2002
  % author: Google
  % ed: J. Leskovec
  % fields: name title A id date author ed kind notes
  % kind: directed graph
  %-------------------------------------------------------------------------------
  ...

and generate it with the Symbol ``:r``::

  julia> matrixdepot("SNAP/web-Google", :r)
  916428x916428 sparse matrix with 5105039 Float64 entries:
	[11343 ,      1]  =  1.0
	[11928 ,      1]  =  1.0
	[15902 ,      1]  =  1.0
	[29547 ,      1]  =  1.0
	[30282 ,      1]  =  1.0
	[31301 ,      1]  =  1.0
	[38717 ,      1]  =  1.0
	[43930 ,      1]  =  1.0
	[46275 ,      1]  =  1.0
	[48193 ,      1]  =  1.0
	[50823 ,      1]  =  1.0
	[56911 ,      1]  =  1.0
	[62930 ,      1]  =  1.0
	[68315 ,      1]  =  1.0
	[71879 ,      1]  =  1.0
	[72433 ,      1]  =  1.0
	[73632 ,      1]  =  1.0
	⋮
	[532967, 916427]  =  1.0
	[547586, 916427]  =  1.0
	[557890, 916427]  =  1.0
	[571471, 916427]  =  1.0
	[580544, 916427]  =  1.0
	[608625, 916427]  =  1.0
	[618730, 916427]  =  1.0
	[622998, 916427]  =  1.0
	[673046, 916427]  =  1.0
	[716616, 916427]  =  1.0
	[720325, 916427]  =  1.0
	[772226, 916427]  =  1.0
	[785097, 916427]  =  1.0
	[788476, 916427]  =  1.0
	[822938, 916427]  =  1.0
	[833616, 916427]  =  1.0
	[417498, 916428]  =  1.0
	[843845, 916428]  =  1.0


Interface to NIST Matrix Market (temporarily suspend)
------------------------------------------------------

Use ``MatrixDepot.get(NAME, collection = :MM)``, where ``NAME`` is
``collection name + '/' + set name + '/' + matrix name`` to download a
test matrix from NIST Matrix Market:
http://math.nist.gov/MatrixMarket/. For example::

  julia> MatrixDepot.get("Harwell-Boeing/lanpro/nos5", collection = :MM)

The way to generate a matrix from NIST Matrix Market and check matrix
information are the same as above.
 

