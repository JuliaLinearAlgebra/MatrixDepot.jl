n = rand(1:7)
A = matrixdepot("invol", n)
# A*A \approx eye(n)
@test_approx_eq_eps A*A eye(n,n) 1e-5
println("'invol' passed test...")
