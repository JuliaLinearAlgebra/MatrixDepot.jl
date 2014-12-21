n = rand(1:10)
A = matrixdepot("randcorr", n)
@test issym(A) # symmetric

# positive semidefinite 
for eig in eigvals(A)
    @test eig >= 0
end

# 1s on the diagonal
@test_approx_eq diag(A) ones(Float64, n)
println("'randcorr' passed test...")
