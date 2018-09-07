n = 8 # rand(2: 10)

# test Hessenberg

A = matrixdepot("rohess", n)

B = tril(A) - diagm(0 => diag(A)) - diagm(-1 => diag(A, -1))

@test B ≈ zeros(n,n) atol=1.0e-7

# test orithogonality

@test A'*A ≈ Matrix(1.0I, n, n) atol=1.0e-7

println("'rohess' passed test...")
