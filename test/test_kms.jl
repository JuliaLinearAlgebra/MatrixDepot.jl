n = 10 # rand(1:10)

A = matrixdepot("kms", n)

B = [1:n;]*ones(n)'
B = abs.(B - B')
B = 0.5.^B

@test A ≈ B

# test complex number

C = [1 2im -4 -8im; -2im 1 2im -4; -4 -2im 1 2im; 8im -4 -2im 1]

@test matrixdepot("kms", 4, 2im) ≈ C

println("'kms' passed test...")
