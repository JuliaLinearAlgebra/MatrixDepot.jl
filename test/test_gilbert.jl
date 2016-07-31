n = rand(1:10)
A = matrixdepot("gilbert", n)
@test issym(A)

B = matrixdepot("gilbert", n, 0.2)
C = matrixdepot("gilbert", Int, n, 0.5)
@test nnz(B) <= nnz(C)
println("'gilbert' passed test ...")
