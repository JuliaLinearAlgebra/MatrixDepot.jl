n = rand(1:10)
@test matrixdepot("cauchy", n) == matrixdepot("cauchy", Float64[1:n;])
x = ones(Float64, n) 
y = x + 2
A = zeros(Float64, n, n)
fill!(A, 0.25)
@test_approx_eq matrixdepot("cauchy", x, y) A
println("'cauchy' passed test...")
