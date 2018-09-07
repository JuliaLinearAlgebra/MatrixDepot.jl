n = 10 # rand(1:10)
A = matrixdepot("smallworld", n)
B = matrixdepot("smallworld", Int, n)

C = matrixdepot("erdrey", n)
C2 = MatrixDepot.shortcuts(C, 0.6)
@test nnz(C) <= nnz(C2)
println("'erdrey' passed test...")
