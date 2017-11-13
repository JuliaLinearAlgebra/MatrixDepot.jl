A = matrixdepot("toeplitz", [1,2,3,4,5])
@test A == matrixdepot("toeplitz", 5)

B = matrixdepot("toeplitz", [1,2,3], [1,4,5])
@test B[1,3] == 5
@test B[3,1] == 3
@test B[1,2] == 4
println("'toeplitz' passed test...")

n = rand(2:10)
A = matrixdepot("prolate", n)
@test issymmetric(A)

println("'prolate' passed test...")
