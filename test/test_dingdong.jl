n = rand(1:10)
p = -2*Float64[1:n;] + (n + 1.5);
C = matrixdepot("cauchy", p)
@test_approx_eq matrixdepot("dingdong", n) C
println("'dingdong' passed test...")
