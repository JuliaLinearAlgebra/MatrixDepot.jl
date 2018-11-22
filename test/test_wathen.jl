n = 5 # rand(1:5)
A = matrixdepot("wathen", n)

@test typeof(A) <: SparseMatrixCSC
@test issymmetric(Matrix(A))
@test isposdef(Matrix(A))

println("'wathen' passed test...")
