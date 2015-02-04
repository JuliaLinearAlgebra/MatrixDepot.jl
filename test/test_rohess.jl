n = rand(2: 10)

# test Hessenberg

A = matrixdepot("rohess", n)

B = tril(A) - diagm(diag(A)) - diagm(diag(A, -1), -1)

@test_approx_eq_eps B zeros(n,n) 1.0e-7

# test orithogonality

@test_approx_eq_eps A'*A eye(n) 1.0e-7

println("'rohess' passed test...")
