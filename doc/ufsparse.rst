.. _ufsparse:

Interface to the UF Sparse Matrix Collection
---------------------------------------------

Before downloading test matrices, we should first update the database::

  julia> updatesparse()

Use ``downloadsparse`` to download a test matrix from the
UF Sparse Matrix Collection: http://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html.
For example::

  julia> downloadsparse("HB/illc1850")

When download is complete, we can generate it using::

  julia> matrixdepot("illc1850")

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
