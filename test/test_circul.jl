n = rand(1:10)
A = ones(Float64, n, n)
@test matrixdepot("circul", ones(Float64,n)) ≈ A

x = [1:n;];
A = matrixdepot("circul", x)
# test symmetry
@test A + A' == (A + A')'
B = matrixdepot("circul", n)
C = matrixdepot("circul", [1:n;], n)
@test B ≈ C

println("'circul' passed test..." )
