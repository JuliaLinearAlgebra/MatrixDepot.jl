n = rand(1: 10)
@test matrixdepot("moler", n) == matrixdepot("moler", n, -1)

M = matrixdepot("moler", n)
for i = 1:n, j = 1:n
    i != j ? (@test_approx_eq M[i,j] min(i,j) - 2) : (@test_approx_eq M[i,i] i)
end

println("'moler' passed test...")
