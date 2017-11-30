A = matrixdepot("hankel", [1:5;])
@test A == matrixdepot("hankel", 5)

r1 = rand(5)
r2 = rand(5)
B = matrixdepot("hankel", r1, r2)
@test issymmetric(B)

C = matrixdepot("hankel", Int, [1,2,3], [3, 4, 5, 6])
@test C[3,4] == 6
@test C[2,4] == C[3,3]
println("'hankel' passed test...")
