n = rand(2:10)

A = matrixdepot("wilkinson", n)

@test A == matrixdepot("wilkinson", Float64, n)

@test issym(full(A))

@test_approx_eq diag(full(A),1) ones(Float64, n - 1)

# symmetric along antidiagonal
@test issym(full(A)[n:-1:1,:])

θ = matrixdepot("wilkinson", 1)

println("'wilkinson' passed test...")


