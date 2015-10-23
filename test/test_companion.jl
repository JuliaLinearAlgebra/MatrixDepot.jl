n = rand(1:10)

A = matrixdepot("companion", Float32, n)
@test diag(A, -1) == ones(Float32, n-1)

v = [3., 4, 2, 10]
B = matrixdepot("companion", v)
@test B[:,end] == v

println("'companion' passed test...")
