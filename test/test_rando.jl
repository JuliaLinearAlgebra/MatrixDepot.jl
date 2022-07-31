let i, n, A, B, C, θ
n = 5 # rand(1:5)

A = matrixdepot("rando", Int, n, 1)
B = matrixdepot("rando", Int, n, 2)
C = matrixdepot("rando", Int, n, 3)

for i in eachindex(A)
    @test A[i] == 0 || A[i] == 1
end

@test all([bi == -1 || bi == 1 for bi in B])
@test all([bi == -1 || bi == 0 || bi == 1 for bi in C])

@test_throws ArgumentError matrixdepot("rando", 3, 3, 5)

θ = matrixdepot("rando", n)
end
println("'rando' passed test...")
