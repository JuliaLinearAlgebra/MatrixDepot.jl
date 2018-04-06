
n = rand(1:10)
p = -2*Float64[1:n;] .+ (n + 1.5);
C = matrixdepot("cauchy", p)
@test matrixdepot("dingdong", n) ≈ C

println("'dingdong' passed test...")
