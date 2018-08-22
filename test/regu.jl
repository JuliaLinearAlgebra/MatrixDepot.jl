n = rand(1:10)

r = matrixdepot("deriv2", n, false)
rf32 = matrixdepot("deriv2", Float32, n, false)

r2 = matrixdepot("deriv2", Float32, n, 2, false)
if mod(n, 2) == 0
    r3 = matrixdepot("deriv2", Float32, n, 3, false)
end
@test_throws ArgumentError matrixdepot("deriv2", Float64, n, 4, false)

A = matrixdepot("deriv2", n)

@test r.A*r.x ≈ r.b
@test issymmetric(r.A)


r = matrixdepot("shaw", 2*n, false)
A = matrixdepot("shaw", 2*n)

print(r)

@test issymmetric(r.A)

rf32 = matrixdepot("shaw", Float32, 2*n, false)

r = matrixdepot("wing", n, false)
A = matrixdepot("wing", n)

@test r.A == matrixdepot("wing", n, 1/3, 2/3, false).A

r = matrixdepot("foxgood", n, false)
rf32 = matrixdepot("foxgood", Float32, n, false)
A = matrixdepot("foxgood", n)

r = matrixdepot("heat", 2*n, false)
rf32 = matrixdepot("heat", Float32, 2*n)

r = matrixdepot("baart", 2*n, false)
rf32 = matrixdepot("baart", Float32, 2*n, false)
A = matrixdepot("baart", 2*n)

n = rand(3:10)
r = matrixdepot("phillips", 4*n, false)
rf32 = matrixdepot("phillips", Float32, 4*n)

A = matrixdepot("gravity", n)

@test matrixdepot("gravity", n, 1, 0, 1, 0.25) == A

r1 = matrixdepot("gravity", Float32, n, 1, false)
r2 = matrixdepot("gravity", Float32, n, 2, false)
r3 = matrixdepot("gravity", Float32, n, 3, false)

A = matrixdepot("blur", n)

@test matrixdepot("blur", n, 3, 0.7) == A

r1 = matrixdepot("blur", Float32, n, false)

A = matrixdepot("spikes", n)

@test matrixdepot("spikes", n, 5) == A

n = rand(5:10)
r1 = matrixdepot("spikes", Float32, n, false)

A = matrixdepot("ursell", n)
@test matrixdepot("ursell", Float64, n) == A
r1 = matrixdepot("ursell", Float32, n, false)
print(r1)

A = matrixdepot("parallax", n)

m, n = size(A)
@test m == 26
@test matrixdepot("parallax", Float64, n) == A
r1 = matrixdepot("parallax", Float32, n, false) 
