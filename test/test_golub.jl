n = rand(6:10)

A = matrixdepot("golub", n)
@test cond(A) > 1.e8

B = matrixdepot("golub", Int, n)
@test cond(A) > 1.e8

println("'golub' passed test...")
