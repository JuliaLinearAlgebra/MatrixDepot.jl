n = rand(1:10)
A = matrixdepot("chebspec", n)
# no null vector for 1-by-1 chebspec.
if n == 1
    n = 2
end
# A has null vector ones(n)
@test_approx_eq_eps A*ones(n) zeros(n) 1e-6
B = matrixdepot("chebspec", n, 1)
# B's eigenvalues have negative real parts.
v = real(eigvals(B))
for i = 1:n
    @test v[i] < 0
end
println("'chebspec' passed test...")
