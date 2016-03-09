.. _interfaces:

Interface to Test Collections
=============================

Before downloading test matrices, it is recommended to first update
the database::

  julia> MatrixDepot.update()


Interface to the UF Sparse Matrix Collection
---------------------------------------------

Use ``matrixdepot(NAME, :get)``, where ``NAME`` is ``collection_name
+'/' + matrix_name``,  to download a test matrix from the
`UF Sparse Matrix Collection <http://www.cise.ufl.edu/research/sparse/matrices/list_by_id.html>`_.
For example::

  julia> matrixdepot("SNAP/web-Google", :get)

.. note:: 
   ``matrixdepot()`` displays all the matrices in the
   collection, including the newly downloaded matrices. All the matrix 
   data can be found by ``matrixdepot("data")``. 
	  
If the matrix name is unique in the collections, we could use
``matrixdepot(matrix_name, :get)`` to download the data. If more than
one matrix has the same name, a list of options will be returned. For
example::
  
  julia> matrixdepot("epb0", :get)
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
  100 83244  100 83244    0     0   109k      0 --:--:-- --:--:-- --:--:--  133k
  download:/home/weijian/.julia/v0.4/MatrixDepot/data/uf/Averous/epb0.tar.gz
  epb0/epb0.mtx

  julia> matrixdepot("1138_bus", :get)
  Try MatrixDepot.get(`name`), where `name` is one of the elements in the following Array:
  2-element Array{AbstractString,1}:
  "HB/1138_bus"                    
  "Harwell-Boeing/psadmit/1138_bus"

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

and generate it with the Symbol ``:r`` or  ``:read`` ::

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

The metadata of a given matrix can be obtained by 
``matrixdepot(collection_name/matrix_name, :read, meta = true)``. For example::

  julia> matrixdepot("TKK/t520", :get)
  julia> matrixdepot("TKK/t520", :read, meta = true)
  Dict{AbstractString,Any} with 3 entries:
    "t520"       => 5563x5563 Symmetric{Float64,SparseMatrixCSC{Float64,Int64}}:…
    "t520_b"     => "%%MatrixMarket matrix array real general\n%-----------------…
    "t520_coord" => "%%MatrixMarket matrix array real general\n%-----------------…

We can use ``matrixdepot(collection_name/*)`` to download all the matrices
in a given ``collection_name``. For example, we can get all the 
matrices contributed by The Mathworks, Inc. by ``matrixdepot("MathWorks/*", :get)``::

   julia> matrixdepot()

   Matrices:
    1) baart            2) binomial         3) blur             4) cauchy        
    5) chebspec         6) chow             7) circul           8) clement       
    9) companion       10) deriv2          11) dingdong        12) fiedler       
    13) forsythe        14) foxgood         15) frank           16) golub         
    17) gravity         18) grcar           19) hadamard        20) hankel        
    21) heat            22) hilb            23) invhilb         24) invol         
    25) kahan           26) kms             27) lehmer          28) lotkin        
    29) magic           30) minij           31) moler           32) neumann       
    33) oscillate       34) parter          35) pascal          36) pei           
    37) phillips        38) poisson         39) prolate         40) randcorr      
    41) rando           42) randsvd         43) rohess          44) rosser        
    45) sampling        46) shaw            47) spikes          48) toeplitz      
    49) tridiag         50) triw            51) ursell          52) vand          
    53) wathen          54) wilkinson       55) wing            56) MathWorks/Harvard500
    57) MathWorks/Kaufhold
    58) MathWorks/Kuu   59) MathWorks/Muu   60) MathWorks/Pd    61) MathWorks/Pd_rhs
    62) MathWorks/pivtol
    63) MathWorks/QRpivot
    64) MathWorks/Sieber
    65) MathWorks/tomography
    66) MathWorks/TS  
   Groups:
    all           data          eigen         ill-cond    
    inverse       pos-def       random        regprob     
    sparse        symmetric  

   julia> matrixdepot("data")
   11-element Array{AbstractString,1}:
   "MathWorks/Harvard500"
   "MathWorks/Kaufhold"  
   "MathWorks/Kuu"       
   "MathWorks/Muu"       
   "MathWorks/Pd"        
   "MathWorks/Pd_rhs"    
   "MathWorks/pivtol"    
   "MathWorks/QRpivot"   
   "MathWorks/Sieber"    
   "MathWorks/tomography"
   "MathWorks/TS" 



Interface to NIST Matrix Market
-------------------------------

Use ``matrixdepot(NAME, :get)``, where ``NAME`` is
``collection name + '/' + set name + '/' + matrix name`` to download a
test matrix from NIST Matrix Market:
http://math.nist.gov/MatrixMarket/. For example::

  julia> matrixdepot("Harwell-Boeing/lanpro/nos5", :get)

If the matrix name is unique, we could also use ``matrixdepot(matrix name, :get)``::

  julia> matrixdepot("bp__1400", :get)
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
  100 28192  100 28192    0     0   4665      0  0:00:06  0:00:06 --:--:-- 10004
  download:/home/weijian/.julia/v0.4/MatrixDepot/data/mm/Harwell-Boeing/smtape/bp__1400.mtx.gz

 
Checking matrix information and generating matrix data are similar to 
the above case::

  julia> matrixdepot("Harwell-Boeing/smtape/bp__1400")

  %%MatrixMarket matrix coordinate real general

  use matrixdepot("Harwell-Boeing/smtape/bp__1400", :read) to read the data

  julia> matrixdepot("Harwell-Boeing/smtape/bp__1400", :read)
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
