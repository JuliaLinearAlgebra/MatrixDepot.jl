n = rand(1:10)
@test matrixdepot("lehmer", Float64, n) == matrixdepot("lehmer", n)

A = ones(n, 1) * [1:n]'
A = A./A'
A = tril(A) + tril(A)' - eye(n)

@test_approx_eq matrixdepot("lehmer", n) A
println("'lehmer' passed test...")

