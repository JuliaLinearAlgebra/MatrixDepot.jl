n = rand(1:10)

A = matrixdepot("wilkinson", n)

@test A == matrixdepot("wilkinson", Float64, n)

@test issym(full(A))

@test_approx_eq diag(full(A),1) ones(Float64, n - 1)

println("'wilkson' passed test...")
