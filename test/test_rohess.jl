n = rand(2: 10)

# test Hessenberg

A = matrixdepot("rohess", n)

B = tril(A) - diagm(diag(A)) - diagm(diag(A, -1), -1)

@test B ≈ zeros(n,n) atol=1.0e-7

# test orithogonality

@test A'*A ≈ eye(n) atol=1.0e-7

println("'rohess' passed test...")
