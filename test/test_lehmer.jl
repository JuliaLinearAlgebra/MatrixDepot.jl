n = 10 # rand(1:10)
@test matrixdepot("lehmer", Float64, n) == matrixdepot("lehmer", n)

A = ones(n, 1) * [1:n;]'
A = A./A'
A = tril(A) + tril(A)' - Matrix(1.0I, n, n)

@test matrixdepot("lehmer", n) ≈ A
println("'lehmer' passed test...")

