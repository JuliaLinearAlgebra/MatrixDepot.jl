n = rand(1:5)
A = matrixdepot("wathen", n)

@test typeof(A) <: SparseMatrixCSC
@test issym(full(A))
@test isposdef(full(A))

println("'wathen' passed test...")
