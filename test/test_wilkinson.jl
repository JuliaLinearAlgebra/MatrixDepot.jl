n = rand(2:10)

A = matrixdepot("wilkinson", n)

@test A == matrixdepot("wilkinson", Float64, n)

@test issymmetric(full(A))

@test diag(full(A),1) ≈ ones(Float64, n - 1)

# symmetric along antidiagonal
@test issymmetric(full(A)[n:-1:1,:])

θ = matrixdepot("wilkinson", 1)

println("'wilkinson' passed test...")


