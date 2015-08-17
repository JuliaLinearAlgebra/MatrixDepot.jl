n = rand(1:10)

r = matrixdepot("deriv2", n)

@test_approx_eq r.A*r.x  r.b
