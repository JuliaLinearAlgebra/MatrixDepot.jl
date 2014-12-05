n = rand(1:10)
A = matrixdepot("pascal", n)

@test isposdef(A)
@test issym(A)

P = diagm((-1).^[0:n-1])
P[:,1] = ones(n, 1)

for j = 2: n-1
    for i = j+ 1:n
        P[i,j] = P[i-1, j] - P[i-1, j-1]
    end
end

@test_approx_eq A P*P'
println("'pascal' passed test...")
