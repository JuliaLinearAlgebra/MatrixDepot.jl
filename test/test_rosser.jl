A = matrixdepot("rosser", Int, 8, 2, 1)

# B is the test matrix used in the paper
#. Rosser, Lanczos, Hestenes and Karush, J. Res. Natl. Bur. Stand. Vol. 47 (1951), pp291-297. 
B = [611 196 -192 407 -8 -52 -49 29;
     196 899 113 -192 -71 -43 -8 -44;
     -192 113 899 196 61 49 8 52;
     407 -192 196 611 8 44 59 -23;
     -8 -71 61 8 411 -599 208 208;
     -52 -43 49 44 -599 411 208 208;
     -49 -8 8 59 208 208 99 -911;
     29 -44 52 -23 208 208 -911 99];

@test A == B

# test eigenvalues
e1 = eigvals(matrixdepot("rosser", 2))
e2 = [500, 510]
@test_approx_eq e1 e2

e1 = eigvals(matrixdepot("rosser", 4))
e2 = [1., 10000, 10199, 10200]
@test_approx_eq_eps e1 e2 1e-2
