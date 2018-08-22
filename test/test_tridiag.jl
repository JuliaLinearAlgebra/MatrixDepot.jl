n = rand(1:10)
@test matrixdepot("tridiag", [1,1,1], [1,1,1,1], [1,1,1]) == matrixdepot("tridiag", 4, 1,1,1)
@test matrixdepot("tridiag", Float64, n) == matrixdepot("tridiag", n)

A = matrixdepot("tridiag", n)
@test isposdef(Array(A))
@test issymmetric(Array(A))
println("'tridiag' passed test...")
