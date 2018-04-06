n = rand(1:7)
A = matrixdepot("invol", n)
# A*A ≈ eye(n)
@test A*A ≈ Matrix(1.0I, n, n) atol = 1e-5
println("'invol' passed test...")
