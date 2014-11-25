# Matrix Depot for Julia

Matrix Depot is a multi-language collection of test matrices.


## Install

To install the package, type

```julia
julia> Pkg.clone("git://github.com/weijianzhang/MatrixDepot.jl.git")
```

## Usage

First load the package:

```julia
julia> using MatrixDepot
```

The only function will be exported is `matrixdepot`.

```julia
julia>? matrixdepot
INFO: Loading help data...
matrixdepot (generic function with 8 methods)
```

Every matrix in the collection is represented by a string `matrix_name`, for
example, the Cauchy matrix is represented by `"cauchy"` and the Hilbert matrix
is represented by `"hilb"`.

The properties of the matrices in the collection are also symbolized by strings
`propertry_name`. For example, the class of the symmetric matrices is symbolized
by `"symmetric"`.

* `matrixdepot(matrix_name, p1, p2, ...)` returns a matrix specified by the
query `matrix_name`. `matrix_name` is a query string. `p1, p2, ...` are input
parameters depending on `matrix_name`.

* `matrixdepot(matrix_name)` returns the parameter options and the properties of
`matrix_name`.

* `matrixdepot(property_name)` returns a list of matrices with the property
`property_name`.

## Examples

We can generate a Hilbert matrix of size 4 by typing


    matrixdepot("hilb", 4)




    4x4 Array{Float64,2}:
     1.0       0.5       0.333333  0.25    
     0.5       0.333333  0.25      0.2     
     0.333333  0.25      0.2       0.166667
     0.25      0.2       0.166667  0.142857



and generate a circul matrix of size 5 by


    matrixdepot("circul", 5)




    5x5 Array{Float64,2}:
     1.0  2.0  3.0  4.0  5.0
     5.0  1.0  2.0  3.0  4.0
     4.0  5.0  1.0  2.0  3.0
     3.0  4.0  5.0  1.0  2.0
     2.0  3.0  4.0  5.0  1.0



If we are not sure about the paramters or properties of a matrix, we can type
the matrix name


    matrixdepot("hilb")

    Hilbert matrix: 
                  
     Input options: 
                  
     (type), dim: the dimension of the matrix
                  
     (type), row_dim, col_dim: the row and column dimension 
                  
     ['inverse', 'ill-cond', 'symmetric', 'pos-def']



    matrixdepot("hadamard")

    Hadamard matrix: 
                  
     Input options: 
                  
     (type), n::Int: the dimension of the matrix, n is a power of 2 
                  
     ['inverse', 'orthogonal', 'eigen']


From the information given by `matrixdepot("hilb")`, we notice we create
a 4-by-6 rectanglular Hilbert matrix by


    matrixdepot("hilb", 4, 6)




    4x6 Array{Float64,2}:
     1.0       0.5       0.333333  0.25      0.2       0.166667
     0.5       0.333333  0.25      0.2       0.166667  0.142857
     0.333333  0.25      0.2       0.166667  0.142857  0.125   
     0.25      0.2       0.166667  0.142857  0.125     0.111111



We can aslo specify the data type


    matrixdepot("hilb", Float16, 5, 3)




    5x3 Array{Float16,2}:
     1.0      0.5      0.33325
     0.5      0.33325  0.25   
     0.33325  0.25     0.19995
     0.25     0.19995  0.16663
     0.19995  0.16663  0.14282



By typing a matrix name, we can see what properties that matrix have.
Conversely, if we type a property, we can see all the matrices (in the
collection) having that property.


    matrixdepot("symmetric")




    3-element Array{ASCIIString,1}:
     "hilb"  
     "cauchy"
     "circul"




    matrixdepot("ill-cond")




    2-element Array{ASCIIString,1}:
     "hilb"  
     "cauchy"




    matrixdepot("inverse")




    3-element Array{ASCIIString,1}:
     "hilb"    
     "hadamard"
     "cauchy"  



Given a property, we can loop through all the matrices having this propery


    # Multiply all matrices of the class "inverse"
    A = eye(4)
    print("Identity matrix")
    for mat in matrixdepot("inverse")
        print(" x $mat matrix")
        A = A * matrixdepot(mat, 4)    
    end
    println(" is equal to")
    A    

    Identity matrix x hilb matrix x hadamard matrix x cauchy matrix is equal to





    4x4 Array{Float64,2}:
     1.54861   1.09306   0.849802  0.696925
     0.833056  0.578056  0.444722  0.362123
     0.578056  0.3975    0.304087  0.246647
     0.444722  0.304087  0.231797  0.187548


