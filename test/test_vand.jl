n = rand(1:10)
m = rand(1:n)
@test matrixdepot("vand", Int, n) == matrixdepot("vand", [1:n;])

A = matrixdepot("vand", Int, n)
@test A[:, m] == [1:n;].^(m-1)
println("'vand' passed test...")
