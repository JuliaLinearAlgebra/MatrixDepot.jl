n = rand(1:10)
A = matrixdepot("erdrey", n)
B = matrixdepot("erdrey", n, ceil(Int, n*log(n)/2))
@test nnz(A) == nnz(B)
@test issym(A)
println("'erdrey' passed test...")

