n = rand(1:10)

r = matrixdepot("deriv2", n)
rf32 = matrixdepot("deriv2", Float32, n)

@test_approx_eq r.A*r.x  r.b
@test issym(r.A)


r = matrixdepot("shaw", 2*n)

@test issym(r.A)

rf32 = matrixdepot("shaw", Float32, 2*n)

r = matrixdepot("wing", n)

@test r.A == matrixdepot("wing", n, 1/3, 2/3).A

r = matrixdepot("foxgood", n)
rf32 = matrixdepot("foxgood", Float32, n)

r = matrixdepot("heat", n)
rf32 = matrixdepot("heat", Float32, n)
