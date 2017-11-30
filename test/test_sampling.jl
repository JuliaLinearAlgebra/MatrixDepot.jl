n = rand(1:7)
@test matrixdepot("sampling", n) == matrixdepot("sampling", Float64, n)

A = matrixdepot("sampling", n)
v = sort(eigvals(A))
@test v ≈ [0: n - 1;] atol=1e-10
println("'sampling' passed test...")
