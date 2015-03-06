.. _ufsparse:

Interface to the UF Sparse Matrix Collection
---------------------------------------------

Before downloading test matrices, we should first update the database::

  julia> updatesparse()

Use ``downloadsparse`` to download a test matrix from the
`UF Sparse Matrix Collection <http://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html>`_.
For example::

  julia> downloadsparse("HB/illc1850")

When download is complete, we can generate it using::

  julia> matrixdepot("illc1850")
  1850x712 sparse matrix with 8636 Float64 entries:
        [1   ,    1]  =  0.27735
	[3   ,    1]  =  0.27735
	[26  ,    1]  =  0.27735
	[28  ,    1]  =  0.27735
	[164 ,    1]  =  0.27735
	[166 ,    1]  =  0.27735
	[1259,    1]  =  0.27735
	[1262,    1]  =  0.27735
	[1277,    1]  =  0.27735
	[1279,    1]  =  0.27735
	[1490,    1]  =  0.27735
	[1686,    1]  =  0.27735
	[1827,    1]  =  0.27735
	[2   ,    2]  =  0.5
	[4   ,    2]  =  0.5
	[5   ,    2]  =  0.5
	[7   ,    2]  =  0.5
	⋮
	[1833,  712]  =  0.0622221
	[1834,  712]  =  0.0626713
	[1835,  712]  =  0.0619168
	[1836,  712]  =  0.0608667
	[1837,  712]  =  0.0610634
	[1838,  712]  =  0.0615807
	[1839,  712]  =  0.0616349
	[1840,  712]  =  0.061873
	[1841,  712]  =  0.061889
	[1842,  712]  =  0.0618741
	[1843,  712]  =  0.0618311
	[1844,  712]  =  0.0618285
	[1845,  712]  =  0.0619223
	[1846,  712]  =  0.061784
	[1847,  712]  =  0.0617665
	[1848,  712]  =  0.0615756
	[1849,  712]  =  0.0616386
	[1850,  712]  =  0.0616394

and check matrix information using::

  julia> matrixdepot("illc1850", :info)
  Dict{ASCIIString,Any} with 10 entries:
  "name"   => "HB/illc1850"
  "A"      => …
  "author" => "M. Saunders"
  "kind"   => "least squares problem"
  "Zeros"  => …
  "b"      => [64.06762598…
  "title"  => "UNSYMMETRIC LEAST-SQUARES PROBLEM.                  SAUNDERS 1979."
  "id"     => 170.0
  "date"   => "1979"
  "ed"     => "I. Duff, R. Grimes, J. Lewis"
