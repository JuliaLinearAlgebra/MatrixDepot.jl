n = rand(1:10)

r = matrixdepot("deriv2", n, false)
rf32 = matrixdepot("deriv2", Float32, n, false)

r2 = matrixdepot("deriv2", Float32, n, 2, false)
if mod(n, 2) == 0
    r3 = matrixdepot("deriv2", Float32, n, 3, false)
end

@test_approx_eq r.A*r.x  r.b
@test issym(r.A)


r = matrixdepot("shaw", 2*n, false)

@test issym(r.A)

rf32 = matrixdepot("shaw", Float32, 2*n, false)

r = matrixdepot("wing", n, false)

@test r.A == matrixdepot("wing", n, 1/3, 2/3, false).A

r = matrixdepot("foxgood", n, false)
rf32 = matrixdepot("foxgood", Float32, n, false)

r = matrixdepot("heat", 2*n)
rf32 = matrixdepot("heat", Float32, 2*n)

r = matrixdepot("baart", 2*n)
rf32 = matrixdepot("baart", Float32, 2*n)

n = rand(1:10)
r = matrixdepot("phillips", 4*n)
rf32 = matrixdepot("phillips", Float32, 4*n)

A = matrixdepot("gravity", n)

@test matrixdepot("gravity", n, 1, 0, 1, 0.25) == A

r1 = matrixdepot("gravity", Float32, n, 1, false)
r2 = matrixdepot("gravity", Float32, n, 2, false)
r3 = matrixdepot("gravity", Float32, n, 3, false)

A = matrixdepot("blur", n)

@test matrixdepot("blur", n, 3, 0.7) == A

r1 = matrixdepot("blur", Float32, n, false)
