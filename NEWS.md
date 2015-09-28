## Matrix Depot Release Notes

v0.3.3 (under preparation)
---------------------------

* display the matrices in group lists alphabetically.

* add group `"all"` for all the matrices in the collection.

v0.3.2 
-------

* make `matrixdepot()` display better information.

* rename `@addproperty` to `@addgroup` and rename `@rmproperty` to `@rmgroup`.

* rename property `regu` to `regprob`.

* set Float64 as the default data type for all the parameterized matrices in
  the collection. 

* update matrix references.
	

v0.3.1 
------

* Add new test problems for regularization methods:

  - `heat`
  - `phillips`
  - `baart`

v0.3.0 
------

* Add test problems for regularization methods:

  - `deriv2`
  - `foxgood`
  - `shaw`
  - `wing`

* Reintroduce Matrix Market interface

* Define new functions for test collection interfaces

  - `matrixdepot(name, :get)`: download matrix data `name`.
  - `matrixdepot(name, :read)`: read matrix data `name`.

v0.2.8
------

* Add new test matrices

  - `oscillate`: a random test matrix for numerical regularization methods.
  - `prolate`: a symmtric ill-conditioned Toeplitz matrix.
  - `toeplitz`: Toeplitz matrix.


v0.2.7
------

* Fix some typos and v0.4 deprecation warnings


v0.2.6
------

* add reference information for test matrices

* update demo


v0.2.5
------

* support accessing matrices by number and UnitRange

* matrices in the collection are ordered alphabetically 

* temporarily suspend the NIST Matrix Market interface

* use base Matrix Market reader for Julia 0.4

v0.2.2
------

* Include an interface to NIST Matrix Market

* reimplement ransvd to include rectangular case

v0.2.1 
------

* Include an interface to the UF Sparse Matrix Collection


v0.1.3
------

* New matrices

    wathen: a sparse symmetric positive random matrix arose from the
    finite element method

* Style the output information


v0.1.2 
------

* New matrices

	rohess: random orthogonal upper Hessenberg matrix
		
	kms: Kac-Murdock-Szego Toeplitz matrix

* Others

	fix test error for randsvd. 

	

v0.1.1 
------

* New matrices 

	wilkinson: Wilkinson's eigenvalue test matrix. 

	rando: random matrix with entries -1, 0 or 1.

	randsvd: random matrix with pre-assigned singular values.

* New property

	random: the matrix is random.

* Optimized Vandermonde matrix generation for better performance

	Thanks to @synapticarbors and @stevengj


v0.1.0
------
* New matrices

	rosser: a matrix with close eigenvalues.

	sampling: a matrix with application in sampling theory.

* Sphinx documentation


TODO
----

* add a subcollection for symmetric quasi-definite linear system?






  
