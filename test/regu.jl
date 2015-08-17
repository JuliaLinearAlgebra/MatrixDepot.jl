n = rand(1:10)

r = matrixdepot("deriv2", n)
rf32 = matrixdepot("deriv2", Float32, n)

@test_approx_eq r.A*r.x  r.b
@test issym(r.A)


r = matrixdepot("shaw", 2*n)

@test issym(r.A)

rf32 = matrixdepot("shaw", Float32, 2*n)
