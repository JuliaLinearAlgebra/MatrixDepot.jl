n = rand(2:10)

kappa = rand(1:6)
mode = rand(1:5)


A = matrixdepot("randsvd", n, kappa, mode)

if mode == 5
    @test cond(A) <= kappa
else 
    @test_approx_eq cond(A) kappa
end

println("'randsvd' passed test...")
