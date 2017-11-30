n = rand(1:10)
m = rand(1:10)
k = rand(1:n)
alpha = 3.0
@test matrixdepot("triw", Float64, m, n, alpha, k) == matrixdepot("triw", m, n, alpha, k)
@test matrixdepot("triw", n) == matrixdepot("triw", Float64, n, n, -1, n-1)

A = matrixdepot("triw", Float32, n)
@test istriu(A)
# B = A'*A has diagnoal B(i,i) = i  
@test diag(A'*A) ≈ [1. : n;]
println("'triw' passed test...")
