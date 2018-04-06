n = rand(1:38)
@test matrixdepot(n) != nothing
@test matrixdepot(1:n) != nothing
