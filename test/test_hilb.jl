n = 10 # rand(1:10)
@test matrixdepot("hilb", n) == matrixdepot("hilb", n, n)
@test matrixdepot("hilb" , Float16, n) ≈ matrixdepot("hilb", Float32, n) 
x = Float64[1:n;]
y = x .- 1
@test matrixdepot("hilb", n) ≈ matrixdepot("cauchy", x, y)
println("'hilb' passed test...")
