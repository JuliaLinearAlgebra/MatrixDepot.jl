n = 8 # rand(1:10)
@test matrixdepot("frank", n) == matrixdepot("frank", Float64, n, 0)

A = zeros(Float64, n, n)
for i = 1:n, j = 1:n
    A[i,j] = i<=j ? n + 1 -j : j == i - 1 ? n - j : zero(Float64)
end
@test matrixdepot("frank", n) ≈ A

A = zeros(Float64, n, n)
for i = 1:n, j = 1:n
    A[n+1-j,n+1-i] = i<=j ? n + 1 -j : j == i - 1 ? n - j : zero(Float64)
end
@test matrixdepot("frank", n, 1) ≈ A

println("'frank' passed test...")
