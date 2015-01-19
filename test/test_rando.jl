n = rand(1:5)

A = matrixdepot("rando", Int, n, 1)
B = matrixdepot("rando", Int, n, 2)
C = matrixdepot("rando", Int, n, 3)

for i in length(A)
    @test A[i] == 0 || A[i] == 1
end

for i in length(B)
    @test B[i] == -1 || B[i] == 1
end

for i in length(C)
    @test C[i] == -1 || C[i] == 0 || C[i] == 1
end

println("'rando' passed test...")
