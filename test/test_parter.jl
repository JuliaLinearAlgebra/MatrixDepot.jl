n = rand(1:10)
@test matrixdepot("parter", n) == matrixdepot("parter", Float64, n)

A = matrixdepot("cauchy", [1.:n;] .+ 0.5, -[1.:n;])

@test matrixdepot("parter", n) ≈ A
println("'parter' passed test...")
