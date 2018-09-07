n = 10 # rand(1:10)
@test matrixdepot("chow", n) == matrixdepot("chow", Float64, n)
@test matrixdepot("chow",Int, 5, 5, 0)[5,1] == 3125
@test diag(matrixdepot("chow", n, 2, 4)) == ones(Float64, n) * 6 
println("'chow' passed test...")
