n = 10 #rand(1:10)
m = n ÷ 2 # rand(1:n)
@test matrixdepot("vand", Int, n) == matrixdepot("vand", [1:n;])

A = matrixdepot("vand", Int, n)
@test A[:, m] == [1:n;].^(m-1)
println("'vand' passed test...")
