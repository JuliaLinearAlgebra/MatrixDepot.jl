n = rand(1:10)
@test matrixdepot("binomial", n) == matrixdepot("binomial", Float64, n)

A = matrixdepot("binomial", n)

@test A*A ≈ 2^(n-1)*Matrix(1.0I, n, n)
println("'binomial' passed test...")
