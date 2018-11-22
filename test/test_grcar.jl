n = 7 # rand(3:10)
@test matrixdepot("grcar", n, 3) == matrixdepot("grcar", n)

A = matrixdepot("grcar", n, n-1)
@test istril(A - ones(n,n))
@test istriu(A + diagm(-1 => ones(n-1)))
println("'grcar' passed test...")
