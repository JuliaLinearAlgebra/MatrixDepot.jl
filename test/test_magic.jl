n = rand(1:10)
M = matrixdepot("magic", n)

@test sum(M,1) == sum(M,1)'
# diagnoal == antidiagnoal
p = [n:-1:1]
@test sum(diag(M)) == sum(diag(M[:,p]))
