n = 9 # rand(1:10)
@test matrixdepot("clement", Float64, n) == matrixdepot("clement", n)

A = matrixdepot("clement", n)
B = matrixdepot("clement", n, 1)

@test diag(A+A', 1) == n*ones(n-1)
@test issymmetric(Array(B))
θ = matrixdepot("clement", 1)
println("'clement' passed test...")
