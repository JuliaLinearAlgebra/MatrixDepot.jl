n = rand(2:10)

kappa = rand(2:6)
mode = rand(1:5)


A = matrixdepot("randsvd", n, kappa, mode)
B = matrixdepot("randsvd", n, kappa)

if mode == 5
    @test cond(A) <= kappa
else 
    @test_approx_eq cond(A) kappa
end

A5 = matrixdepot("randsvd", n, kappa, 5)
A4 = matrixdepot("randsvd", n, kappa, 4)
A3 = matrixdepot("randsvd", n, kappa, 3)
A2 = matrixdepot("randsvd", n, kappa, 2)
A1 = matrixdepot("randsvd", n, kappa, 1)
θ = matrixdepot("randsvd", 1)
try 
    matridepot("randsvd", 5, kappa, 6)
catch ArgumentError
    println("randsvd: invalid mode value")
end

println("'randsvd' passed test...")
