n = 10 # rand(1:10)
m = 9  # rand(1:10)
@test matrixdepot("kahan", n) == matrixdepot("kahan", n, 1.2, 25.)
@test matrixdepot("kahan", m ,n, 4, 5) == matrixdepot("kahan", Float64, m, n, 4, 5)

@test matrixdepot("kahan", n, 0.5*pi, 1) ≈ Matrix(1.0I, n,n) 
println("'kahan' passed test...")
