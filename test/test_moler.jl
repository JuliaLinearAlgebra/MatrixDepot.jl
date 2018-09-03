n = 10 # rand(1: 10)
@test matrixdepot("moler", n) == matrixdepot("moler", n, -1)

M = matrixdepot("moler", n)
for i = 1:n, j = 1:n
    if i != j
        @test M[i,j] ≈ min(i,j) - 2
    else
        @test M[i,i] ≈ i
    end
end

println("'moler' passed test...")
