n = rand(40:100)
A = matrixdepot("gilbert", n)
@test issymmetric(A)

B = matrixdepot("gilbert", n, 0.2)
C = matrixdepot("gilbert", Int, n, 0.5)
@test nnz(B) <= nnz(C)
println("'gilbert' passed test ...")
