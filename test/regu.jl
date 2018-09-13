n = 7 # rand(1:10)

r = mdopen("deriv2", n, false)
@test r !== nothing
rf32 = mdopen("deriv2", Float32, n, false)
@test rf32 !== nothing

r2 = mdopen("deriv2", Float32, n, 2, false)
@test r2 !== nothing
@test mdopen("deriv2", Float32, (n ÷ 2) * 2, 3, false) !== nothing

@test_throws ArgumentError matrixdepot("deriv2", Float64, n, 4, false)

A = matrixdepot("deriv2", n)
@test A == r.A

@test r.A*r.x ≈ r.b
@test issymmetric(r.A)


r = mdopen("shaw", 2*n, false)
@test r !== nothing
A = matrixdepot("shaw", 2*n)
@test A !== nothing

# print(r)

@test issymmetric(r.A)

rf32 = matrixdepot("shaw", Float32, 2*n, false)
@test rf32 !== nothing

r = mdopen("wing", n, false)
@test r !== nothing
A = matrixdepot("wing", n)
@test A !== nothing && A == r.A

@test r.A == matrixdepot("wing", n, 1/3, 2/3, false)

r = mdopen("foxgood", n, false)
@test r !== nothing
rf32 = mdopen("foxgood", Float32, n, false)
@test rf32 !== nothing
A = matrixdepot("foxgood", n)
@test A !== nothing

r = mdopen("heat", 2*n, false)
@test r !== nothing
rf32 = matrixdepot("heat", Float32, 2*n)
@test rf32 !== nothing

r = mdopen("baart", 2*n, false)
@test r !== nothing
rf32 = mdopen("baart", Float32, 2*n, false)
@test rf32 !== nothing
A = matrixdepot("baart", 2*n)
@test A !== nothing

n = 8 # rand(3:10)
r = mdopen("phillips", 4*n, false)
@test r !== nothing
rf32 = matrixdepot("phillips", Float32, 4*n)
@test rf32 !== nothing

A = matrixdepot("gravity", n)
@test A !== nothing

@test matrixdepot("gravity", n, 1, 0, 1, 0.25) == A

@test mdopen("gravity", Float32, n, 1, false) !== nothing
@test mdopen("gravity", Float32, n, 2, false) !== nothing
@test mdopen("gravity", Float32, n, 3, false) !== nothing
@test_throws ErrorException matrixdepot("gravity", Float32, n, 4, false)

A = matrixdepot("blur", n)
@test A !== nothing

@test matrixdepot("blur", n, 3, 0.7) == A

r1 = mdopen("blur", Float32, n, false)
@test r1 !== nothing

A = matrixdepot("spikes", n)
@test A !== nothing

@test matrixdepot("spikes", n, 5) == A

n = 7 # rand(5:10)
r1 = mdopen("spikes", Float32, n, false)
@test r1 !== nothing

A = matrixdepot("ursell", n)
@test A !== nothing
@test matrixdepot("ursell", Float64, n) == A
r1 = mdopen("ursell", Float32, n, false)
@test r1 !== nothing
# print(r1)

A = matrixdepot("parallax", n)
@test A !== nothing

m, n = size(A)
@test m == 26
@test matrixdepot("parallax", Float64, n) == A
r1 = mdopen("parallax", Float32, n, false) 
@test r1 !== nothing
