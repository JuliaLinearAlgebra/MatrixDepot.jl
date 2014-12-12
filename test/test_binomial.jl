n = rand(1:10)
@test matrixdepot("binomial", n) == matrixdepot("binomial", Int, n)

A = matrixdepot("binomial", n)

@test_approx_eq A*A 2^(n-1)*eye(n)
println("'binomial' passed test...")
