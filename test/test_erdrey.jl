n = rand(1:10)
A = matrixdepot("erdrey", n)
B = matrixdepot("erdrey", n, ceil(Int, n*log(n)/2))
@test nnz(A) == nnz(B)
@test issym(A)
m = rand(1:n)
C = matrixdepot("erdrey", n, m)
@test nnz(C) == m*2
println("'erdrey' passed test...")
