n = 10 # rand(1:10)
A = matrixdepot("erdrey", n)
B = matrixdepot("erdrey", n, ceil(Int, n*log(n)/2))
@test nnz(A) == nnz(B)
@test issymmetric(A)

C = matrixdepot("erdrey", 10, 3)
@test nnz(C) == 6
println("'erdrey' passed test...")
