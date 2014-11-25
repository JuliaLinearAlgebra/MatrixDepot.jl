for n in [2, 4, 16]
    H = matrixdepot("hadamard", Int, n)
    @test H*H' == n*eye(n)
end
println("'hadamard' passed test...")
