
for n in [1, 3, 4, 5, 6]
    M = matrixdepot("magic", n)
    @test sum(M, dims=1) == sum(M, dims=2)'
    k = rand(1:n)
    @test sum(M, dims=1)[k] == sum(M, dims=2)'[k]
    # diagonal == antidiagnoal
    p = [n:-1:1;]
    @test sum(diag(M)) == sum(diag(M[:,p]))
end
@test size(matrixdepot("magic", 2)) == (2, 2)

println("'magic' passed test...")
