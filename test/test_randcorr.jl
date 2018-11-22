n = 10 # rand(1:10)
A = matrixdepot("randcorr", n)
@test issymmetric(A) # symmetric

# positive semidefinite 
for eig in eigvals(A)
    @test eig >= 0
end

# 1s on the diagonal
@test diag(A) ≈ ones(Float64, n)
println("'randcorr' passed test...")
