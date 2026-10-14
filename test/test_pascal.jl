n = 10 # rand(1:10)
A = matrixdepot("pascal", n)

@test isposdef(A)
@test issymmetric(A)

P = diagm(0 => (-1).^[0:n-1;])
P[:,1] = ones(n, 1)

for j = 2: n-1
    for i = j+ 1:n
        P[i,j] = P[i-1, j] - P[i-1, j-1]
    end
end

@test A ≈ P*P'
@test matrixdepot("pascal", 36)[36, 36] ≈ binomial(BigInt(70), 35) # Int overflow
println("'pascal' passed test...")
