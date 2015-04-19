n = rand(1:10)
@test matrixdepot("fiedler", n) == matrixdepot("fiedler", [1:n;])

A = matrixdepot("fiedler", n)

@test issym(A)
for i = 1:n, j =1:n
    @test_approx_eq A[i,j] abs(i-j)
end
println("'fiedler' passed test...")
