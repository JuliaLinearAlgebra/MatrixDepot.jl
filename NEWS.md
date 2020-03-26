## Matrix Depot Release Notes

v1.0.0

* additional meta-data for suite-sparse collection specially SVD data
* improved testing
* ready for windows
* adapt documentation


v0.8.1
--------

* fix problem with serialization of database
* adapt data recognition for SuiteSparse index file
* fix typo in sample problem "baart"

v0.8.0
--------

* Improved API (see: README.md and `?MatrixDepot`)
* Keyword search
* Logical expressions of predicates
* Access current remote repositories (2018.8)
* some bugs fixed - see git log
* Documentation adapted

V0.7.0
--------

* Drop support for Julia v0.6

v0.6.0
--------

* Adapt to Juliav0.7 - deprecations and minor bugs

v0.5.6
--------

* Fix various typos in documentation, thanks to @jiahao.

v0.5.5
--------

* Improve the documentation of regularization problems.

v0.5.4
--------

* Fix `SparseMatrix` deprecation for Julia v0.5.

* Enhance the UF sparse matrix collection's meta data storage.

v0.5.3
---------

* Add new test problems for regularization methods:

 - `parallax`: http://matrixdepotjl.readthedocs.org/en/latest/regu.html#term-parallax

* Allow the UF sparse matrix collection interface to store meta data. 

v0.5.2
------

* Update documentation.

v0.5.1
------

* Fix a bug in user defined generators.

* Fix reference formatting.

* Allow accessing matrices by different subtypes of Integer or UnitRange.

v0.5.0 
------

* Drop support for Julia 0.3

* Various enhancements, including:

  - better formatted documentation for matrix generators.

  - automatically includes all Julia files in the directory
	`myMatrixDepot`, so that adding new matrix generators is more
	convenient and flexible. See
	http://matrixdepotjl.readthedocs.org/en/latest/user.html more
	details.

v0.4.3
------

* Ignore any changes to the directory `myMatrixDepot`. This means
  any local changes to `myMatrixDepot\generator.jl` won't be affected when
  we upgrade to a new version of Matrix Depot.

* Add more tests

v0.4.2
------

* Fix: build failing on Julia v0.5

* Prevent `@rmgroup` from introducing extra empty lines

* Add more tests

v0.4.1
------

* Add new feature:

   - allow users to download a group of matrices from UF sparse matrix collection.

v0.4.0 
------

* Add new feature:

	- allow users to add new matrix generators.
	  see http://matrixdepotjl.readthedocs.org/en/latest/user.html

* Add new test problems for regularization methods:

	- `spikes`
	- `ursell`

* Add new test matrices:

	- `golub` see http://blogs.mathworks.com/cleve/2015/08/24/golub-matrices-deceptively-ill-conditioned/
    - `companion`: Companion matrix.

* Document function `matrixdepot`.

v0.3.5 
----------

* Add new test problems for regularization methods:

	- `gravity`
	- `blur`

* Add new feature:

	- access matrices by a mixture of number and range. For example,
	  we could do `matrixdepot(1:4, 6, 10:15)`. 

* Implement all three examples for `deriv2`.

v0.3.4 
-------

* fix Julia v0.5 String deprecation.

* implement a more general matrix generator. 

* add new matrix: Hankel matrix `hankel` 

* enumerate matrix input options

* add option `matrixonly` for regularization problems.

v0.3.3 
-------

* display the matrices in group lists alphabetically.

* add group `"all"` for all the matrices in the collection.

* fix pascal matrix overflow error.

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
  - `prolate`: a symmetric ill-conditioned Toeplitz matrix.
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






  
