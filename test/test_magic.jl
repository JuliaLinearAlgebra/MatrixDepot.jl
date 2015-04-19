n = rand(1:10)
M = matrixdepot("magic", n)
if n == 2 # no magic square of order 2
    return
else
    @test sum(M,1) == sum(M,2)'
    k = rand(1:n)
    @test sum(M,1)[k] == sum(M,2)'[k]
    # diagnoal == antidiagnoal
    p = [n:-1:1;]
    @test sum(diag(M)) == sum(diag(M[:,p]))
end

println("'magic' passed test...")
