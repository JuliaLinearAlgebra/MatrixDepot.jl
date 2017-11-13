n = rand(1:10)
@test matrixdepot("forsythe", n) == 
                matrixdepot("forsythe", n, sqrt(eps(Float64)), zero(Float64))
@test matrixdepot("forsythe", n , 1, 2) == matrixdepot("forsythe", Float64, n, 1, 2)
T = Float32
a = convert(T, 3)
b = convert(T, 4)
@test matrixdepot("forsythe", T, n, 3, 4) == matrixdepot("forsythe", n, a, b)

if n != 1 # this test does not apply for n = 1
    A = zeros(T, n, n)
    for i = 1:n, j = 1:n
        A[i,j] = j == i ? b : 
                  j == i + 1 ? 1.0 : 
                  (i == n && j == 1)? a : 0.0
    end

    @test matrixdepot("forsythe", n, a, b) ≈ A
end
println("'forsythe' passed test...")

