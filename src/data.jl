matrixdict = @compat Dict("hilb" => hilb, "hadamard" => hadamard,
                          "cauchy" => cauchy, "circul" => circul,
                          "dingdong" => dingdong, "frank" => frank,
                          "invhilb" => invhilb, "forsythe" => forsythe,
                          "magic" => magic, "grcar" => grcar,
                          "triw" => triw, "moler" => moler,
                          "pascal" => pascal, "kahan" => kahan,
                          "pei" => pei, "vand" => vand,
                          "invol" => invol, "chebspec" => chebspec,
                          "lotkin" => lotkin, "clement" => clement,
                          "fiedler" => fiedler, "minij" => minij,
                          "binomial" => binomialm, "tridiag" => tridiag,
                          "lehmer" => lehmer, "parter" => parter,
                          "chow" => chow, "randcorr" => randcorr,
                          "poisson" => poisson, "neumann" => neumann,
                          "rosser" => rosser, "sampling" => sampling,
                          "wilkinson" => wilkinson, "rando" => rando,
                          "randsvd" => randsvd, "rohess" => rohess,
                          "kms" => kms, "wathen" => wathen,
                          "oscillate" => oscillate, "toeplitz" => toeplitz,
                          "prolate" => prolate, "deriv2" => deriv2,
                          "shaw" => shaw, "wing" => wing,
                          "foxgood" => foxgood, "heat" => heat,
                          "baart" => baart, "phillips" => phillips
                          );

matrixinfo =
@compat Dict("hilb" => "Hilbert matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix
             \n [type,] row_dim, col_dim: the row and column dimension
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference: M.-D. Choi, Tricks or treats with the Hilbert matrix,
             Amer. Math. Monthly, 90 (1983), pp. 301-312.
             N.J. Higham, Accuracy and Stability of Numerical Algorithms,
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 1996; sec. 26.1.",

             "invhilb" => "Inverse of Hilbert matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix
             \n ['inverse', 'ill-cond', 'symmetric','pos-def']
             \n Reference: M.-D. Choi, Tricks or treats with the Hilbert matrix,
             Amer. Math. Monthly, 90 (1983), pp. 301-312.
             N.J. Higham, Accuracy and Stability of Numerical Algorithms,
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 1996; sec. 26.1.",

             "hadamard" => "Hadamard matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix, n is a power of 2
             \n ['inverse', 'orthogonal', 'eigen']
             \n Reference: S.W. Golomb and L.D. Baumert, The search for
             Hadamard matrices, Amer. Math. Monthly, 70 (1963) pp. 12-17.",

             "cauchy" => "Cauchy matrix:
             \n Input options:
             \n vec1, vec2: two vectors
             \n vec: a vector
             \n [type,] dim: the dimension of the matrix
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference:  N.J. Higham, Accuracy and Stability of Numerical Algorithms,
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 1996; sec. 26.1.",

             "circul" => "Circul matrix:
             \n Input options:
             \n vec, n: a vector and the column dimension
             \n vec: a vector
             \n [type,] dim: the dimension of the matrix
             \n ['symmetric', 'pos-def', 'eigen']
             \n Reference:  P.J. Davis, Circulant Matrices, John Wiley, 1977.",

             "dingdong" => "Dingdong matrix:
             \n Input options:
             \n [type,] n: the dimension of the matrix.
             \n ['symmetric', 'eigen']
             \n Reference: J.C. Nash, Compact Numerical Methods for
             Computers: Linear Algebra and Function Minimisation,
             second edition, Adam Hilger, Bristol, 1990 (Appendix 1).",

             "frank" => "Frank matrix:
             \n Input options:
             \n [type,] n, k: n is the dimension of the matrix, k = 0 or 1.
             If k = 1 the matrix reflect about the anti-diagonal.
             \n [type,] n: n is the dimension of the matrix.
             \n ['ill-cond', 'eigen']
             \n Reference:  W.L. Frank, Computing eigenvalues of complex matrices
             by determinant evaluation and by methods of Danilewski and Wielandt,
             J. Soc. Indust. Appl. Math., 6 (1958), pp. 378-392 (see pp. 385, 388).",

             "forsythe" => "Forsythe matrix:
             \n Input options:
             \n [type,] n, alpha, lambda: n is the dimension of the matrix.
             alpha and lambda are scalars.
             \n [type,] n: alpha = sqrt(eps(type)) and lambda = 0.
             \n ['inverse', 'ill-cond', 'eigen']
             \n This generator is adapted from Nicholas J. Higham's Test Matrix Toolbox.",

             "magic" => "Magic square matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix.
             \n ['inverse']",

             "grcar" => "Grcar Matrix:
             \n Input options:
             \n [type,] dim, k: dim is the dimension of the matrix and
             k is the number of superdiagonals.
             \n [type,] dim: the dimension of the matrix.
             \n ['eigen']
             \n Reference: J.F. Grcar, Operator coefficient methods
             for linear equations, Report SAND89-8691, Sandia National
             Laboratories, Albuquerque, New Mexico, 1989 (Appendix 2).",

             "triw" => "Triw Matrix:
             \n Input options:
             \n [type,] row_dim, col_dim, alpha, k: row_dim and col_dim
             are row and column dimension of the matrix. alpha is a
             scalar representing the entries on the superdiagonals.
             k is the number superdiagonals.
             \n [type,] dim
             \n ['inverse', 'ill-cond']
             \n Reference:  G.H. Golub and J.H. Wilkinson, Ill-conditioned
             eigensystems and the computation of the Jordan canonical form,
             SIAM Review, 18(4), 1976, pp. 578-6",

             "moler" => "Moler Matrix:
             \n Input options:
             \n [type,] dim, alpha: dim is the dimension of the matrix,
             alpha is a scalar.
             \n [type,] dim: alpha = -1.
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference:  J.C. Nash, Compact Numerical Methods for Computers:
             Linear Algebra and Function Minimisation, second edition,
             Adam Hilger, Bristol, 1990 (Appendix 1).",

             "pascal" => "Pascal Matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix.
             \n ['Inverse', 'ill-cond', 'symmetric', 'pos-def', 'eigen']
             \n Reference: R. Brawer and M. Pirovino, The linear algebra of
             the Pascal matrix, Linear Algebra and Appl., 174 (1992),
             pp. 13-23 (this paper gives a factorization of L = PASCAL(N,1)
             and a formula for the elements of L^k).
             N.J. Higham, Accuracy and Stability of Numerical Algorithms
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 1996; sec. 26.4.",

             "kahan" => "Kahan Matrix:
             \n Input options:
             \n [type,] m, n, theta, pert: m, n are the row and column
             dimensions of the matrix. theta and pert are scalars.
             \n [type,] dim, theta, pert: dim is the dimension of the matrix.
             \n [type,] dim: theta = 1.2, pert = 25.
             \n ['inverse', 'ill-cond']
             \n Reference: W. Kahan, Numerical linear algebra, Canadian Math.
             Bulletin, 9 (1966), pp. 757-801." ,

             "pei" => "Pei Matrix:
             \n Input options:
             \n [type,] dim, alpha: dim is the dimension of the matrix.
             alpha is a scalar.
             \n [type,] dim
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference: M.L. Pei, A test matrix for inversion procedures,
             Comm. ACM, 5 (1962), p. 508.",

             "vand" => "Vandermonde Matrix:
             \n Input options:
             \n vec, dim: vec is a vector, dim is the number of columns.
             \n vec
             \n [type,] dim
             \n ['inverse', 'ill-cond']
             \n Reference: N.J. Higham, Stability analysis of algorithms
             for solving confluent Vandermonde-like systems, SIAM J.
             Matrix Anal. Appl., 11 (1990), pp. 23-41.",

             "invol" => "Involutory Matrix:
             \n Input options:
             \n [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'ill-cond', 'eigen']
             \n Reference: A.S. Householder and J.A. Carpenter, The
             singular values of involutory and of idempotent matrices,
             Numer. Math. 5 (1963), pp. 234-237.",

             "chebspec" => "Chebyshev spectral differentiation matrix:
             \n Input options:
             \n [type,] dim, k: dim is the dimension of the matrix and
             k = 0 or 1.
             \n [type,] dim
             \n ['eigen']
             \n Reference: L.N. Trefethen and M.R. Trummer, An instability
             phenomenon in spectral methods, SIAM J. Numer. Anal., 24 (1987), pp. 1008-1023.",

             "lotkin" => "Lotkin Matrix:
             \n Input options:
             \n [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'ill-cond', 'eigen']
             \n Reference: M. Lotkin, A set of test matrices, MTAC, 9 (1955), pp. 153-161.",

             "clement" => "Clement Matrix:
             \n Input options:
             \n [type,] dim, k: dim is the dimension of the matrix.
             If k = 0, the matrix is Tridiagonal. If k = 1, the matrix
             is SymTridiagonal.
             \n [type,] dim: k = 0
             \n ['inverse', 'symmetric', 'eigen']
             \n Reference: P.A. Clement, A class of triple-diagonal
             matrices for test purposes, SIAM Review, 1 (1959), pp. 50-52.",

             "fiedler" => "Fiedler Matrix:
             \n Input options:
             \n vec: vec is a vector.
             \n [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'symmetric', 'eigen']
             \n Reference: G. Szego, Solution to problem 3705, Amer. Math.
              Monthly, 43 (1936), pp. 246-259.
             J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
             Birkhauser, Basel, and Academic Press, New York, 1977, p. 159.",

             "minij" => "MIN[I,J] Matrix:
             \n Input options:
             \n [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'symmetric', 'pos-def', 'eigen']
             \n Reference: J. Fortiana and C. M. Cuadras, A family of matrices,
             the discretized Brownian bridge, and distance-based regression,
             Linear Algebra Appl., 264 (1997), 173-188.  (For the eigensystem of A.)",

             "binomial" => "Binomial Matrix:
             \n Input options:
             \n [type,] dim: dim is the dimension of the matrix.
             \n ['']",

             "tridiag" => "Tridiagonal Matrix:
             \n Input options:
             \n v1, v2, v3: v1 and v3 are sub- superdiagonal vectors.
             v2 is the diagonal vector.
             \n [type,] dim, x, y, z: dim is the dimension of the matrix. x, y, z are
             scalars. x and z are sub- superdiagonal elements, y is diagonal
             element.
             \n [type,] dim: x = -1, y = 2, z = -1.
             \n ['inverse', 'ill-cond', 'pos-def', 'eigen']
             \n Reference: J. Todd, Basic Numerical Mathematics, Vol. 2:
             Numerical Algebra, Birkhauser, Basel, and Academic Press,
             New York, 1977, p. 155.",

             "lehmer" => "Lehmer Matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix.
             \n ['inverse', 'symmetric', 'pos-def']
             \n Reference: M. Newman and J. Todd, The evaluation of
             matrix inversion programs, J. Soc. Indust. Appl. Math.,
             6 (1958), pp. 466-476.
             Solutions to problem E710 (proposed by D.H. Lehmer): The inverse
             of a matrix, Amer. Math. Monthly, 53 (1946), pp. 534-535.",

             "parter" => "Parter Matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix.
             \n ['eigen']
             \n Reference: The MathWorks Newsletter, Volume 1, Issue 1,
             March 1986, page 2. S.V. Parter, On the distribution of the
             singular values of Toeplitz matrices, Linear Algebra and
             Appl., 80 (1986), pp. 115-130.",

             "chow" => "Chow Matrix:
             \n Input options:
             \n [type,] dim, alpha, delta: dim is dimension of the matrix.
             alpha, delta are scalars such that A[i,i] = alpha + delta and
             A[i,j] = alpha^(i + 1 -j) for j + 1 <= i.
             \n [type,] dim: alpha = 1, delta = 0.
             \n ['eigen']
             \n Reference: T.S. Chow, A class of Hessenberg matrices with known
             eigenvalues and inverses, SIAM Review, 11 (1969), pp. 391-395.",

             "randcorr" => "Random Correlation Matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix.
             \n ['symmetric', 'pos-semidef', 'random']",

             "poisson" => "Poisson Matrix:
             \n Input options:
             \n [type,] n: the dimension of the matirx is n^2.
             \n ['inverse', 'symmetric', 'pos-def', 'eigen', 'sparse']
             \n Reference: G.H. Golub and C.F. Van Loan, Matrix Computations,
             second edition, Johns Hopkins University Press, Baltimore,
             Maryland, 1989 (Section 4.5.4).",

             "neumann" => "Neumann Matrix:
             \n [type,] n: the dimension of the matrix is n^2.
             \n ['eigen', 'sparse']
             \n Reference: R.J. Plemmons, Regular splittings and the
             discrete Neumann problem, Numer. Math., 25 (1976), pp. 153-161.",

             "rosser" => "Rosser Matrix:
             \n Input options:
             \n [type,] dim, a, b: dim is the dimension of the matrix.
             dim must be a power of 2.
             a and b are scalars. For dim = 8, a = 2 and b = 1, the generated
             matrix is the test matrix used by Rosser.
             \n [type,] dim: a = b = rand(1:5)
             \n ['eigen', 'ill-cond', 'random']
             \n Reference: J.B. Rosser, C. Lanczos, M.R. Hestenes, W. Karush,
             Separation of close eigenvalues of a real symmetric matrix,
             Journal of Research of the National Bureau of Standards, v(47)
             (1951)",

             "sampling" => "Matrix with Application in Sampling Theory:
             \n Input options:
             \n vec: vec is a vector with no repeated elements.
             \n [type,] dim: the dimension of the matrix.
             \n ['eigen']
             \n Reference: L. Bondesson and I. Traat, A nonsymmetric matrix
             with integer eigenvalues, linear and multilinear algebra, 55(3)
             (2007), pp. 239-247",

             "wilkinson" => "Wilkinson Matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix.
             \n ['symmetric', 'eigen']
             \n Reference: J.H. Wilkinson, Error analysis of direct methods
             of matrix inversion, J. Assoc. Comput. Mach., 8 (1961),  pp. 281-330.",

             "rando" => "Random Matrix with Element -1, 0, 1:
             \n Input options:
             \n [type,] m, n, k: m and n are row and column dimensions,
             k = 1: entries are 0 or 1.
             k = 2: entries are -1 or 1.
             k = 3: entries are -1, 0 or 1.
             \n [type,] n, k
             \n [type,] n: k = 1
             \n ['random']",

             "randsvd" => "Random Matrix with Pre-assigned Singular Values:
             \n Input options:
             \n [type,] m, n, kappa, mode: m, n are the dimensions of the matrix.
             kappa is the condition number of the matrix.
             mode = 1: one large singular value.
             mode = 2: one small singular value.
             mode = 3: geometrically distributed singular values.
             mode = 4: arithmetrically distributed singular values.
             mode = 5: random singular values with  unif. dist. logarithm.
             \n [type,] n, kappa, mode: m = n
             \n [type,] n, kappa: mode = 3
             \n [type,] n: kappa = sqrt(1/eps()), mode = 3.
             \n ['ill-cond', 'random']
             \n Reference: N.J. Higham, Accuracy and Stability of Numerical
             Algorithms, Society for Industrial and Applied Mathematics,
             Philadelphia, PA, USA, 1996; sec. 26.3.",

             "rohess" => "Random Orthogonal Upper Hessenberg Matrix:
             \n Input options:
             \n [type,] n : n is the dimension of the matrix.
             \n ['random']
             \n Reference:  W.B. Gragg, The QR algorithm for unitary
             Hessenberg matrices, J. Comp. Appl. Math., 16 (1986), pp. 1-8.",

             "kms" => "Kac-Murdock-Szego Toeplitz Matrix:
             \n Input:
             \n [type,] n, rho: n is the dimension of the matrix, rho is a
             scalar such that A[i,j] = rho^(abs(i-j)).
             \n [type,] n: rho = 0.5
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference: W.F. Trench, Numerical solution of the eigenvalue
             problem for Hermitian Toeplitz matrices, SIAM J. Matrix Analysis
             and Appl., 10 (1989), pp. 135-146 (and see the references therein).",

             "wathen" => "Wathen Matrix:
             \n Input options:
             \n [type,] nx, ny: the dimension of the matrix is equal to
             3 * nx * ny + 2 * nx * ny + 1;
             \n [type,] n: nx = ny = n.
             \n ['symmetric', 'pos-def', 'eigen', 'random', 'sparse']
             \n Reference: A.J. Wathen, Realistic eigenvalue bounds for
             the Galerkin mass matrix, IMA J. Numer. Anal., 7 (1987),
             pp. 449-457.",

             "oscillate" => "Oscillating Matrix:
             \n A matrix A is called oscillating if A is totally 
             nonnegative and if there exists an integer q > 0 such that 
             A^q is totally positive.
             \n Input options:
             \n Σ: the singular vaule spectrum of the matrix;
             \n [type,] n, mode: n is the dimension of the matrix. 
                 mode = 1: geometrically distributed singular values.
                 mode = 2: arithmetrically distributed singular values.
             \n [type,] n: mode = 1.
             \n ['symmetric','pos-def', 'random', 'eigen'] 
             \n Reference: Per Christian Hansen, Test matrices for 
             regularization methods. SIAM J. SCI. COMPUT Vol 16, 
             No2, pp 506-512 (1995).",

             "toeplitz" => "Toeplitz Matrix:
             \n A Toeplitz matrix is a matrix in which each descending 
             diagonal from left to right is constant.
             \n Input options:
             \n vc, vr: vc and vr are the first row and column of the matrix;
             \n v: symmatric case, i.e., vc = vr = v;
             \n [type,] n: the dimension of the matrix is n, v = [1:n] is the first 
                row and column vector.",
             
             "prolate" => "Prolate Matrix:
             \n A prolate matrix is a symmetirc, ill-conditioned Toeplitz matrix.
             \n Input options:
             \n [type,] n, w: the dimension of the matrix is n, w is a real scalar;
             \n [type,] n: the case when w = 0.25.
             \n Reference: J.M. Varah. The Prolate Matrix. Linear Algebra and Appl.
             187:267--278, 1993.",
             
             "deriv2" => "Computation of the Second Derivative:
             \n A classical test problem for regularization algorithms.
             \n Input options:
             \n [type,] n: the dimension of the matrix is n.
             \n Reference: P.C. Hansen, Regularization tools: A MATLAB pacakge for 
             analysis and solution of discrete ill-posed problems. 
             Numerical Algorithms, 6(1994), pp.1-35",

             "shaw" => "One-Dimensional Image Restoration Model:
             \n This test problem uses a first-kind Fredholm integral equation
             to model a one-dimentional image restoration situation.
             \n Input options:
             \n [type,] n: the dimesion of the matrix n must be even.
             \n Reference: C.B. Shaw, Jr., Improvements of the resolution of 
             an instrument by numerical solution of an integral equation. 
             J. Math. Ana. Appl. 37 (1972), 83-112.",

             "wing" => "A Problem with a Discontinuous Solution:
             \n Input options:
             \n [type,] n, t1, t2: the dimension of matrix is n. t1 and t2 are 
             two real scalars such that 0 < t1 < t2 < 1;
             \n [type,] n: t1 = 1/3 and t2 = 2/3.
             \n Reference: G.M. Wing, A Primer on Integral Equations of the 
             First Kind, SIAM, 1991, p. 109.",

             "foxgood" => "Severely Ill-posed Problem Suggested by Fox & Goodwin:
             This is a model problem discretized by simple quadrature, which does 
             not satifiy the discrete Picard condition for the small singular 
             values.
             \n Input options:
             \n [type,] n: the dimension of the matrix is n.
             \n Reference: C.T. H. Baker, The Numerical Treatment of Integral
             Equations, Clarendon Press, Oxford, 1977, p. 665.",

             "heat" => "Inverse Heat Equation:
             \n Input options:
             \n [type,] n, κ: the dimension of the matrix is n and n 
             must be even. κ controls the ill-conditioning of the matrix. 
             (κ = 5 gives a well-conditioned problem and κ = 1 
              gives an ill conditoned problem);
             \n [type,] n: κ = 1.
             \n Reference: A.S. Carasso, Determining surface temperatures 
             from interior observations, SIAM J. Appl. Math. 42 (1982), 558-574.",

             "baart" => "Fredholm Integral Equation of the Fisrt Kind:
             \n Input options:
             \n [type,] n: the dimenstion of the matrix is n.
             \n Reference: M.L. Baart, The use of auto-correlation for 
             pesudo-rank determination in noisy ill-conditioned linear-squares
             problems, IMA, J. Numer. Anal. 2 (1982), 241-247.",
             
             "phillips" => "Phillips's \"famous\" problem:
             \n Input options:
             \n [type,] n: the dimenstion of the matrix is n.
             \n Reference: D.L. Phillips, A technique for the numerical 
             solution of certain integral equations of the first kind, J. ACM
             9 (1962), 84-97."
             );

matrixclass =
@compat Dict("symmetric" => ["hilb", "cauchy", "circul", "dingdong",
                             "invhilb", "moler", "pascal", "pei",
                             "clement", "fiedler", "minij", "tridiag",
                             "lehmer", "randcorr", "poisson", "wilkinson",
                              "kms", "wathen", "oscillate", "prolate", 
                             "deriv2", "shaw"],
             "inverse" => ["hilb", "hadamard", "cauchy", "invhilb",
                           "forsythe", "magic", "triw", "moler", "pascal",
                           "kahan", "pei", "vand", "invol", "lotkin",
                           "clement", "fiedler", "minij", "tridiag",
                           "lehmer", "poisson", "kms" ],
             "ill-cond" => ["hilb", "cauchy", "frank", "invhilb",
                            "forsythe", "triw", "moler", "pascal",
                            "kahan","pei", "vand", "invol", "lotkin",
                            "tridiag", "rosser", "randsvd", "kms", 
                            "oscillate", "prolate", "deriv2", "shaw", 
                            "wing", "foxgood", "phillips"],
             "pos-def" => ["hilb", "cauchy", "circul", "invhilb",
                           "moler", "pascal", "pei", "minij", "tridiag",
                           "lehmer", "poisson", "kms", "wathen", "oscillate"],
             "eigen" =>   ["hadamard", "circul", "dingdong", "frank",
                           "forsythe", "grcar", "pascal", "invol","chebspec",
                           "lotkin", "clement", "fiedler", "minij",
                           "tridiag", "parter", "chow", "poisson", "neumann",
                           "rosser", "sampling", "wilkinson","wathen", 
                           "oscillate"],

             "sparse" => ["poisson", "neumann", "wathen"],
             "random" => ["rosser", "rando", "randcorr", "randsvd", "rohess",
                          "wathen", "oscillate"],
             "regu" => ["deriv2", "shaw", "wing", "foxgood", "heat", "baart",
                        "phillips"]
               );
