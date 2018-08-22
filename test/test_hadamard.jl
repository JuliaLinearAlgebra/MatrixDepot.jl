let n
for n in [2, 4, 16]
    H = matrixdepot("hadamard", Int, n)
    @test H*H' == n * Matrix(1.0I, n, n)
end

@test_throws ArgumentError matrixdepot("hadamard", 0)
end    
println("'hadamard' passed test...")
