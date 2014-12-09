n = rand(1:10)
m = rand(1:10)
@test matrixdepot("kahan", n) == matrixdepot("kahan", n, 1.2, 25.)
@test matrixdepot("kahan", m ,n, 4, 5) == matrixdepot("kahan", Float64, m, n, 4, 5)

@test_approx_eq matrixdepot("kahan", n, 0.5*pi, 1) eye(n,n) 
println("'kahan' passed test...")
