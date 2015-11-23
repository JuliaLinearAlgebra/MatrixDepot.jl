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
                          "hankel" => hankel, "golub" => golub,
                          "companion" => companion,
                          "prolate" => prolate, "deriv2" => deriv2,
                          "shaw" => shaw, "wing" => wing,
                          "foxgood" => foxgood, "heat" => heat,
                          "baart" => baart, "phillips" => phillips,
                          "gravity" => gravity, "blur" => blur,
                          "spikes" => spikes, "ursell" => ursell
                          );

matrixinfo =
@compat Dict("hilb" => "Hilbert matrix:
             The Hilbert matrix is a very ill conditioned matrix. 
             It is symmetric positive definite and totally positive. 
             \n Input options:
             \n 1. [type,] dim: the dimension of the matrix;
             \n 2. [type,] row_dim, col_dim: the row and column dimension.
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference: M. D. Choi, Tricks or treats with the Hilbert matrix,
             Amer. Math. Monthly, 90 (1983), pp. 301-312.
             N. J. Higham, Accuracy and Stability of Numerical Algorithms,
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 2002; sec. 28.1.",

             "invhilb" => "Inverse of Hilbert matrix:
             \n Input options:
             \n 1. [type,] dim: the dimension of the matrix.
             \n ['inverse', 'ill-cond', 'symmetric','pos-def']
             \n Reference: M. D. Choi, Tricks or treats with the Hilbert matrix,
             Amer. Math. Monthly, 90 (1983), pp. 301-312.
             N. J. Higham, Accuracy and Stability of Numerical Algorithms,
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 2002; sec. 28.1.",

             "companion" => "Companion matrix:
             The companion matrix to a monic polynomial 
                 a(x) = a_0 + a_1x + ... + a_{n-1}x^{n-1} + x^n
             is the n-by-n matrix with ones on the subdiagonal and 
             the last column given by the coefficients of a(x).     
             \n Input options:
             \n 1. [type,] vec: vec is a vector of coefficients;
             \n 2. [type,] n: vec = [1:n;], the dimension of the matrix 
             is n.",

             "hadamard" => "Hadamard matrix:
             The Hadamard matrix is a square matrix whose entries are 
             1 or -1. It was named after Jacques Hadamard. The rows of 
             a Hadamard matrix are orthogonal.
             \n Input options:
             \n 1. [type,] dim: the dimension of the matrix, n is a power of 2.
             \n ['inverse', 'orthogonal', 'eigen']
             \n Reference: S. W. Golomb and L. D. Baumert, The search for
             Hadamard matrices, Amer. Math. Monthly, 70 (1963) pp. 12-17.",

             "cauchy" => "Cauchy matrix:
             Given two vectors x and y, the (i,j) entry of the Cauchy 
             matrix is 1/(x[i]-y[j]). 
             \n Input options:
             \n 1. [type,] vec1, vec2: two vectors;
             \n 2. [type,] vec: a vector;
             \n 3. [type,] dim: the dimension of the matrix.
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference:  N. J. Higham, Accuracy and Stability of Numerical Algorithms,
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 2002; sec. 28.1.",

             "circul" => "Circulant matrix:
             A circulant matrix has the property that each row is obtained
             by cyclically permuting the entries of the previous row one
             step forward.
             \n Input options:
             \n 1. [type,] vec, n: a vector and the column dimension;
             \n 2. [type,] vec: a vector;
             \n 3. [type,] dim: the dimension of the matrix.
             \n ['symmetric', 'pos-def', 'eigen']
             \n Reference:  P. J. Davis, Circulant Matrices, John Wiley, 1977.",

             "dingdong" => "Dingdong matrix:
             The Dingdong matrix is a symmetric Hankel matrix invented
             by DR. F. N. Ris of IBM, Thomas J Watson Research Centre.
             The eigenvalues cluster around π/2 and -π/2.
             \n Input options:
             \n 1. [type,] n: the dimension of the matrix.
             \n ['symmetric', 'eigen']
             \n Reference: J. C. Nash, Compact Numerical Methods for
             Computers: Linear Algebra and Function Minimisation,
             second edition, Adam Hilger, Bristol, 1990 (Appendix 1).",

             "frank" => "Frank matrix:
             The Frank matrix is an upper Hessenberg matrix with 
             determinant 1. The eigenvalues are real, positive and 
             very ill conditioned.
             \n Input options:
             \n 1. [type,] n, k: n is the dimension of the matrix, k = 0 or 1.
             If k = 1 the matrix reflect about the anti-diagonal;
             \n 2. [type,] n: n is the dimension of the matrix.
             \n ['ill-cond', 'eigen']
             \n Reference:  W. L. Frank, Computing eigenvalues of complex matrices
             by determinant evaluation and by methods of Danilewski and Wielandt,
             J. Soc. Indust. Appl. Math., 6 (1958), pp. 378-392 (see pp. 385, 388).",

             "forsythe" => "Forsythe matrix:
             The Forsythe matrix is a n-by-n perturbed Jordan block.
             \n Input options:
             \n 1. [type,] n, alpha, lambda: n is the dimension of the matrix.
             alpha and lambda are scalars;
             \n 2. [type,] n: alpha = sqrt(eps(type)) and lambda = 0.
             \n ['inverse', 'ill-cond', 'eigen']
             \n This generator is adapted from N. J. Higham's Test Matrix Toolbox.",

             "magic" => "Magic square matrix:
             The magic matrix is a matrix with integer entries such that 
             the row elements, column elements, diagonal elements and 
             anti-diagonal elements all add up to the same number.
             \n Input options:
             \n 1. [type,] dim: the dimension of the matrix.
             \n ['inverse']",

             "grcar" => "Grcar Matrix:
             The Grcar matrix is a Toeplitz matrix with sensitive 
             eigenvalues.
             \n Input options:
             \n 1. [type,] dim, k: dim is the dimension of the matrix and
             k is the number of superdiagonals;
             \n 2. [type,] dim: the dimension of the matrix.
             \n ['eigen']
             \n Reference: J. F. Grcar, Operator coefficient methods
             for linear equations, Report SAND89-8691, Sandia National
             Laboratories, Albuquerque, New Mexico, 1989 (Appendix 2).",

             "triw" => "Triw Matrix:
             Upper triangular matrices discussed by Wilkinson and others.
             \n Input options:
              \n 1. [type,] row_dim, col_dim, α, k: row_dim and col_dim
             are row and column dimension of the matrix. α is a
             scalar representing the entries on the superdiagonals.
             k is the number superdiagonals;
             \n 2. [type,] dim: the dimension of the matrix.
             \n ['inverse', 'ill-cond']
             \n Reference:  G. H. Golub and J. H. Wilkinson, Ill-conditioned
             eigensystems and the computation of the Jordan canonical form,
             SIAM Review, 18(4), 1976, pp. 578-6",

             "moler" => "Moler Matrix:
             The Moler matrix is a symmetric positive definite matrix. 
             It has one small eigenvalue.
             \n Input options:
             \n 1. [type,] dim, alpha: dim is the dimension of the matrix,
             alpha is a scalar;
             \n 2. [type,] dim: alpha = -1.
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference:  J.C. Nash, Compact Numerical Methods for Computers:
             Linear Algebra and Function Minimisation, second edition,
             Adam Hilger, Bristol, 1990 (Appendix 1).",

             "pascal" => "Pascal Matrix:
             The Pascal matrix’s anti-diagonals form the Pascal’s triangle.
             \n Input options:
             \n 1. [type,] dim: the dimension of the matrix.
             \n ['Inverse', 'ill-cond', 'symmetric', 'pos-def', 'eigen']
             \n Reference: R. Brawer and M. Pirovino, The linear algebra of
             the Pascal matrix, Linear Algebra and Appl., 174 (1992),
             pp. 13-23 (this paper gives a factorization of L = PASCAL(N,1)
             and a formula for the elements of L^k).
             N. J. Higham, Accuracy and Stability of Numerical Algorithms
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 2002; sec. 28.4.",

             "kahan" => "Kahan Matrix:
             The Kahan matrix is an upper trapezoidal matrix, i.e., the 
             (i,j) element is equal to 0 if i > j. The useful range of 
             θ is 0 < θ < π. The diagonal is perturbed by 
             pert*eps()*diagm([n:-1:1;]).
             \n Input options:
             \n 1. [type,] m, n, θ, pert: m, n are the row and column
             dimensions of the matrix. theta and pert are scalars;
             \n 2. [type,] dim, θ, pert: dim is the dimension of the matrix;
             \n 3. [type,] dim: θ = 1.2, pert = 25.
             \n ['inverse', 'ill-cond']
             \n Reference: W. Kahan, Numerical linear algebra, Canadian Math.
             Bulletin, 9 (1966), pp. 757-801." ,

             "pei" => "Pei Matrix:
             The Pei matrix is a symmetric matrix with known inversion.
             \n Input options:
             \n 1. [type,] dim, alpha: dim is the dimension of the matrix.
             alpha is a scalar;
             \n 2. [type,] dim: the dimension of the matrix.
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference: M. L. Pei, A test matrix for inversion procedures,
             Comm. ACM, 5 (1962), p. 508.",

             "vand" => "Vandermonde Matrix:
             The inverse and determinat are known explicity.
             \n Input options:
             \n 1. [type,] vec, dim: vec is a vector, dim is the number of columns.
             \n 2. [type,] vec
             \n 3. [type,] dim
             \n ['inverse', 'ill-cond']
             \n Reference: N. J. Higham, Stability analysis of algorithms
             for solving confluent Vandermonde-like systems, SIAM J.
             Matrix Anal. Appl., 11 (1990), pp. 23-41.",

             "invol" => "Involutory Matrix:
             An involutory matrix is a matrix that is its own inverse.
             \n Input options:
             \n 1. [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'ill-cond', 'eigen']
             \n Reference: A. S. Householder and J. A. Carpenter, The
             singular values of involutory and of idempotent matrices,
             Numer. Math. 5 (1963), pp. 234-237.",

             "chebspec" => "Chebyshev spectral differentiation matrix:
             If k = 0,the generated matrix is nilpotent and a vector with 
             all one entries is a null vector. If k = 1, the generated 
             matrix is nonsingular and well-conditioned. Its eigenvalues 
             have negative real parts.
             \n Input options:
             \n 1. [type,] dim, k: dim is the dimension of the matrix and
             k = 0 or 1.
             \n 2. [type,] dim
             \n ['eigen']
             \n Reference: L. N. Trefethen and M. R. Trummer, An instability
             phenomenon in spectral methods, SIAM J. Numer. Anal., 24 (1987), pp. 1008-1023.",

             "lotkin" => "Lotkin Matrix:
             The Lotkin matrix is the Hilbert matrix with its first row 
             altered to all ones. It is unsymmetric, ill-conditioned and 
             has many negative eigenvalues of small magnitude.
             \n Input options:
             \n 1. [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'ill-cond', 'eigen']
             \n Reference: M. Lotkin, A set of test matrices, MTAC, 9 (1955), pp. 153-161.",

             "clement" => "Clement Matrix:
             The Clement matrix is a tridiagonal matrix with zero
             diagonal entries. If k = 1, the matrix is symmetric.
             \n Input options:
             \n 1. [type,] dim, k: dim is the dimension of the matrix.
             If k = 0, the matrix is Tridiagonal. If k = 1, the matrix
             is SymTridiagonal;
             \n 2. [type,] dim: k = 0.
             \n ['inverse', 'symmetric', 'eigen']
             \n Reference: P. A. Clement, A class of triple-diagonal
             matrices for test purposes, SIAM Review, 1 (1959), pp. 50-52.",

             "fiedler" => "Fiedler Matrix:
             The Fiedler matrix is symmetric matrix with a dominant 
             positive eigenvalue and all the other eigenvalues are negative.
             \n Input options:
             \n 1. [type,] vec: vec is a vector;
             \n 2. [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'symmetric', 'eigen']
             \n Reference: G. Szego, Solution to problem 3705, Amer. Math.
              Monthly, 43 (1936), pp. 246-259.
             J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
             Birkhauser, Basel, and Academic Press, New York, 1977, p. 159.",

             "minij" => "MIN[I,J] Matrix:
             A matrix with (i,j) entry min(i,j). It is a symmetric positive
             definite matrix. The eigenvalues and eigenvectors are known 
             explicitly. Its inverse is tridiagonal.
             \n Input options:
             \n 1. [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'symmetric', 'pos-def', 'eigen']
             \n Reference: J. Fortiana and C. M. Cuadras, A family of matrices,
             the discretized Brownian bridge, and distance-based regression,
             Linear Algebra Appl., 264 (1997), 173-188.  (For the eigensystem of A.)",

             "binomial" => "Binomial Matrix:
             The matrix is a multiple of an involutory matrix.
             \n Input options:
             \n 1. [type,] dim: dim is the dimension of the matrix.
             \n ['']",

             "tridiag" => "Tridiagonal Matrix:
             \n Input options:
             \n 1. [type,] v1, v2, v3: v1 and v3 are vectors of subdiagonal
             and superdiagonal elements, respectively, and v2 is a vector
             of diagonal elements;
             \n 2. [type,] dim, x, y, z: dim is the dimension of the matrix,
             x, y, z are scalars. x and z are the subdiagonal and 
             superdiagonal elements, respectively, and y is the diagonal 
             elements;
             \n 3. [type,] dim: x = -1, y = 2, z = -1. This matrix is also
             known as the second difference matrix.
             \n ['inverse', 'ill-cond', 'pos-def', 'eigen']
             \n Reference: J. Todd, Basic Numerical Mathematics, Vol. 2:
             Numerical Algebra, Birkhauser, Basel, and Academic Press,
             New York, 1977, p. 155.",

             "lehmer" => "Lehmer Matrix:
             The Lehmer matrix is a symmetric positive definite matrix. 
             It is totally nonnegative. The inverse is tridiagonal and 
             explicitly known
             \n Input options:
             \n 1. [type,] dim: the dimension of the matrix.
             \n ['inverse', 'symmetric', 'pos-def']
             \n Reference: M. Newman and J. Todd, The evaluation of
             matrix inversion programs, J. Soc. Indust. Appl. Math.,
             6 (1958), pp. 466-476.
             Solutions to problem E710 (proposed by D.H. Lehmer): The inverse
             of a matrix, Amer. Math. Monthly, 53 (1946), pp. 534-535.",

             "parter" => "Parter Matrix:
             The Parter matrix is a Toeplitz and Cauchy matrix 
             with singular values near π.
             \n Input options:
             \n 1. [type,] dim: the dimension of the matrix.
             \n ['eigen']
             \n Reference: The MathWorks Newsletter, Volume 1, Issue 1,
             March 1986, page 2. S. V. Parter, On the distribution of the
             singular values of Toeplitz matrices, Linear Algebra and
             Appl., 80 (1986), pp. 115-130.",

             "chow" => "Chow Matrix:
             The Chow matrix is a singular Toeplitz lower Hessenberg matrix.
             \n Input options:
             \n 1. [type,] dim, alpha, delta: dim is dimension of the matrix.
             alpha, delta are scalars such that A[i,i] = alpha + delta and
             A[i,j] = alpha^(i + 1 -j) for j + 1 <= i;
             \n 2. dim: alpha = 1, delta = 0.
             \n ['eigen']
             \n Reference: T. S. Chow, A class of Hessenberg matrices with known
             eigenvalues and inverses, SIAM Review, 11 (1969), pp. 391-395.",

             "randcorr" => "Random Correlation Matrix:
             A random correlation matrix is a symmetric positive 
             semidefinite matrix with 1s on the diagonal.
             \n Input options:
             \n 1. [type,] dim: the dimension of the matrix.
             \n ['symmetric', 'pos-semidef', 'random']",

             "poisson" => "Poisson Matrix:
             A block tridiagonal matrix from Poisson’s equation. 
             This matrix is sparse, symmetric positive definite and 
             has known eigenvalues.
             \n Input options:
             \n 1. [type,] n: the dimension of the matirx is n^2.
             \n ['inverse', 'symmetric', 'pos-def', 'eigen', 'sparse']
             \n Reference: G. H. Golub and C. F. Van Loan, Matrix Computations,
             second edition, Johns Hopkins University Press, Baltimore,
             Maryland, 1989 (Section 4.5.4).",

             "neumann" => "Neumann Matrix:
             A singular matrix from the discrete Neumann problem. 
             This matrix is sparse and the null space is formed by 
             a vector of ones
             \n Input options:
             \n 1. [type,] n: the dimension of the matrix is n^2.
             \n ['eigen', 'sparse']
             \n Reference: R. J. Plemmons, Regular splittings and the
             discrete Neumann problem, Numer. Math., 25 (1976), pp. 153-161.",

             "rosser" => "Rosser Matrix:
             The Rosser matrix’s eigenvalues are very close together 
             so it is a challenging matrix for many eigenvalue algorithms.
             \n Input options:
             \n 1. [type,] dim, a, b: dim is the dimension of the matrix.
             dim must be a power of 2.
             a and b are scalars. For dim = 8, a = 2 and b = 1, the generated
             matrix is the test matrix used by Rosser.
             \n 2. [type,] dim: a = rand(1:5), b = rand(1:5).
             \n ['eigen', 'ill-cond', 'random']
             \n Reference: J. B. Rosser, C. Lanczos, M. R. Hestenes, W. Karush,
             Separation of close eigenvalues of a real symmetric matrix,
             Journal of Research of the National Bureau of Standards, v(47)
             (1951)",

             "sampling" => "Matrix with Application in Sampling Theory:
             A nonsymmetric matrix with eigenvalues 0, 1, 2, ... n-1.
             \n Input options:
             \n 1. [type,] vec: vec is a vector with no repeated elements;
             \n 2. [type,] n: the dimension of the matrix is n.
             \n ['eigen']
             \n Reference: L. Bondesson and I. Traat, A nonsymmetric matrix
             with integer eigenvalues, linear and multilinear algebra, 55(3)
             (2007), pp. 239-247",

             "wilkinson" => "Wilkinson Matrix:
             The Wilkinson matrix is a symmetric tridiagonal matrix with 
             pairs of nearly equal eigenvalues. The most frequently used                 case is matrixdepot(\"wilkinson\", 21).
             \n Input options:
             \n 1. [type,] dim: the dimension of the matrix.
             \n ['symmetric', 'eigen']
             \n Reference: J. H. Wilkinson, Error analysis of direct methods
             of matrix inversion, J. Assoc. Comput. Mach., 8 (1961),  pp. 281-330.",

             "rando" => "Random Matrix with Element -1, 0, 1:
             \n Input options:
             \n 1. [type,] m, n, k: m and n are row and column dimensions,
             k = 1: entries are 0 or 1.
             k = 2: entries are -1 or 1.
             k = 3: entries are -1, 0 or 1.
             \n 2. [type,] n, k: m = n;
             \n 3. [type,] n: k = 1.
             \n ['random']",

             "randsvd" => "Random Matrix with Pre-assigned Singular Values:
             \n Input options:
             \n 1. [type,] m, n, kappa, mode: m, n are the dimensions of the matrix.
             kappa is the condition number of the matrix.
             mode = 1: one large singular value.
             mode = 2: one small singular value.
             mode = 3: geometrically distributed singular values.
             mode = 4: arithmetrically distributed singular values.
             mode = 5: random singular values with  unif. dist. logarithm;
             \n 2. [type,] n, kappa, mode: m = n;
             \n 3. [type,] n, kappa: mode = 3;
             \n 4. [type,] n: kappa = sqrt(1/eps()), mode = 3.
             \n ['ill-cond', 'random']
             \n Reference: N. J. Higham, Accuracy and Stability of Numerical
             Algorithms, Society for Industrial and Applied Mathematics,
             Philadelphia, PA, USA, 2002; sec. 28.3.",

             "rohess" => "Random Orthogonal Upper Hessenberg Matrix:
             The matrix is constructed via a product of Givens rotations.
             \n Input options:
             \n 1. [type,] n : n is the dimension of the matrix.
             \n ['random']
             \n Reference:  W. B. Gragg, The QR algorithm for unitary
             Hessenberg matrices, J. Comp. Appl. Math., 16 (1986), pp. 1-8.",

             "kms" => "Kac-Murdock-Szego Toeplitz Matrix:
             \n Input:
             \n 1. [type,] n, rho: n is the dimension of the matrix, rho is a
             scalar such that A[i,j] = rho^(abs(i-j)).
             \n 2. [type,] n: rho = 0.5
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             \n Reference: W. F. Trench, Numerical solution of the eigenvalue
             problem for Hermitian Toeplitz matrices, SIAM J. Matrix Analysis
             and Appl., 10 (1989), pp. 135-146 (and see the references therein).",

             "wathen" => "Wathen Matrix:
             Wathen Matrix is a sparse, symmetric positive, random matrix 
             arose from the finite element method. The generated matrix is 
             the consistent mass matrix for a regular nx-by-ny grid of 
             8-nodes.
             \n Input options:
             \n 1. [type,] nx, ny: the dimension of the matrix is equal to
             3 * nx * ny + 2 * nx * ny + 1;
             \n 2. [type,] n: nx = ny = n.
             \n ['symmetric', 'pos-def', 'eigen', 'random', 'sparse']
             \n Reference: A. J. Wathen, Realistic eigenvalue bounds for
             the Galerkin mass matrix, IMA J. Numer. Anal., 7 (1987),
             pp. 449-457.",

             "oscillate" => "Oscillating Matrix:
             \n A matrix A is called oscillating if A is totally 
             nonnegative and if there exists an integer q > 0 such that 
             A^q is totally positive.
             \n Input options:
             \n 1. [type,] Σ: the singular vaule spectrum of the matrix;
             \n 2. [type,] n, mode: n is the dimension of the matrix. 
                 mode = 1: geometrically distributed singular values.
                 mode = 2: arithmetrically distributed singular values.
             \n 3. [type,] n: mode = 1.
             \n ['symmetric','pos-def', 'random', 'eigen'] 
             \n Reference: Per Christian Hansen, Test matrices for 
             regularization methods. SIAM J. SCI. COMPUT Vol 16, 
             No2, pp 506-512 (1995).",

             "toeplitz" => "Toeplitz Matrix:
             \n A Toeplitz matrix is a matrix in which each descending 
             diagonal from left to right is constant.
             \n Input options:
             \n 1. [type,] vc, vr: vc and vr are the first column and row of the matrix;
             \n 2. [type,] v: symmatric case, i.e., vc = vr = v;
             \n 3. [type,] n: the dimension of the matrix is n, v = [1:n;] is the first
                row and column vector.",

             "hankel" => "Hankel Matrix:
             \n A Hankel matrix is a matrix that is symmetric and constant 
             across the anti-diagonals.
             \n Input options:
             \n 1. [type,] vc, vr: vc and vc are the first column and last row of the 
             matrix. If the last element of vc differs from the first element 
             of vr, the last element of rc prevails;
             \n 2. [type,] v: vc = vr = v;
             \n 3. [type,] n: the dimension of the matrix is n, v = [1:n;].",
             
             "prolate" => "Prolate Matrix:
             \n A prolate matrix is a symmetirc, ill-conditioned Toeplitz 
             matrix.
             \n Input options:
             \n 1. [type,] n, w: the dimension of the matrix is n, w is a real scalar;
             \n 2. [type,] n: the case when w = 0.25.
             \n Reference: J. M. Varah. The Prolate Matrix. Linear Algebra and Appl.
             187:267--278, 1993.",

             "golub" => "Golub Matrix:
             Golub matrix is the product of two random unit lower and upper
             triangular matrices respectively. LU factorization without pivoting
             fails to reveal that such matrices are badly conditioned.
             \n Input options:
             \n 1. [type,] n: the dimension of the matrix is n.
             \n Reference: D. Viswanath and N. Trefethen. Condition Numbers of 
             Random Triangular Matrices, SIAM J. Matrix Anal. Appl. 19, 564-581,
             1998.",
             
             "deriv2" => "Computation of the Second Derivative:
             \n A classical test problem for regularization algorithms:
             This is a mildly ill-posed problem. It is a discretization
             of a first kind Fredholm integral equation whose kernel K
             is the Green's function for the second derivative. 
             \n Input options:
             \n 1. [type,] n, example, [matrixonly]: the dimension of the 
             matrix is n.  One choose between between the following right-hand
             g and solution f: 
             (a) example = 1 gives g(s) = (s^3 - s)/6, f(t) = t.
             (b) example = 2 gives g(s) = exp(s) + (1 -e)s - 1, f(t) = exp(t)
             (c) example = 3 gives g(s) = | (4s^3 - 3s)/24,  s < 0.5
                                          | (-4s^3 + 12s^2 - 9s + 1)/24, s>= 0.5
                                   f(t) = | t, t < 0.5
                                   g(t) = | 1- t, t >= 0.5.
             If matrixonly = false, the linear system A, b, x will be 
             generated.  (matrixonly = true by default.)             
             \n 2. [type,] n, [matrixonly]: example = 1.
             \n Reference: P. C. Hansen, Regularization tools: A MATLAB pacakge for 
             analysis and solution of discrete ill-posed problems. 
             Numerical Algorithms, 6(1994), pp.1-35",

             "shaw" => "One-Dimensional Image Restoration Model:
             \n This test problem uses a first-kind Fredholm integral equation
             to model a one-dimentional image restoration situation.
             \n Input options:
             \n 1. [type,] n, [matrixonly]: the dimesion of the matrix n must be even.
             If matrixonly = false, the linear system A, b, x will be generated. 
             (matrixonly = true by default.)
             \n Reference: C. B. Shaw, Jr., Improvements of the resolution of 
             an instrument by numerical solution of an integral equation. 
             J. Math. Ana. Appl. 37 (1972), 83-112.",

             "wing" => "A Problem with a Discontinuous Solution:
             \n Input options:
             \n 1. [type,] n, t1, t2, [matrixonly]: the dimension of matrix is n. 
             t1 and t2 are two real scalars such that 0 < t1 < t2 < 1. 
             If matrixonly = false, the linear system A, b, x will be generated. 
             (matrixonly = true by default.);
             \n 2. [type,] n, [matrixonly]: t1 = 1/3 and t2 = 2/3.
             \n Reference: G. M. Wing, A Primer on Integral Equations of the 
             First Kind, SIAM, 1991, p. 109.",

             "foxgood" => "Severely Ill-posed Problem Suggested by Fox & Goodwin:
             This is a model problem discretized by simple quadrature, which does 
             not satifiy the discrete Picard condition for the small singular 
             values.
             \n Input options:
             \n 1. [type,] n, [matrixonly]: the dimension of the matrix is n.
             If matrixonly = false, the linear system A, b, x will be generated. 
             (matrixonly = true by default.)
             \n Reference: C. T. H. Baker, The Numerical Treatment of Integral
             Equations, Clarendon Press, Oxford, 1977, p. 665.",

             "heat" => "Inverse Heat Equation:
             \n Input options:
             \n 1. [type,] n, κ, [matrixonly]: the dimension of the matrix is n and n 
             must be even. κ controls the ill-conditioning of the matrix.
             (κ = 5 gives a well-conditioned problem and κ = 1 
             gives an ill conditoned problem). 
             If matrixonly = false, the linear system A, b, x will be generated. 
             (matrixonly = true by default.);
             \n 2. [type,] n, [matrixonly]: κ = 1.
             \n Reference: A. S. Carasso, Determining surface temperatures 
             from interior observations, SIAM J. Appl. Math. 42 (1982), 558-574.",

             "baart" => "Fredholm Integral Equation of the Fisrt Kind:
             \n Input options:
             \n 1. [type,] n, [matrixonly]: the dimenstion of the matrix is n.
             If matrixonly = false, the linear system A, b, x will be generated. 
             (matrixonly = true by default.)
             \n Reference: M. L. Baart, The use of auto-correlation for 
             pesudo-rank determination in noisy ill-conditioned linear-squares
             problems, IMA, J. Numer. Anal. 2 (1982), 241-247.",
             
             "phillips" => "Phillips's \"famous\" problem:
             \n Input options:
             \n 1. [type,] n, [matrixonly]: the dimenstion of the matrix is n.
             If matrixonly = false, the linear system A, b, x will be generated. 
             (matrixonly = true by default.)
             \n Reference: D. L. Phillips, A technique for the numerical 
             solution of certain integral equations of the first kind, J. ACM
             9 (1962), 84-97.",
             
             "gravity" => "One-dimensional gravity surverying problem:
             Discretization of a 1-D model problem in gravity surveying, in 
             which a mass distribution f(t) is located at depth d, while the
             vertical component of the gravity field g(s) is measured at the
             surface. 
             \n Input options:
             \n 1. [type,] n, example, a, b, d, [matrixonly]: n is the dimension
             of the matrix. Three examples are implemented.
             (a) example = 1 gives f(t) = sin(pi*t) + 0.5*sin(2*pi*t).
             (b) example = 2 gives f(t) = piecewise linear function.
             (c) example = 3 gives f(t) = piecewise constant function. 
             The t integration interval is fixed to [0, 1], while the s 
             integration interval [a, b] can be specified by the user.
             The parameter d is the depth at which the magnetic deposit is 
             located. The larger the d, the faster the decay of the singular 
             values. 
             If matrixonly = false, the linear system A, b, x will be generated.
             (matrixonly = true by default.)
             \n 2. [type,] n, example, [matrixonly]: a = 0, b = 1, d = 0.25;
             \n 3. [type,] n, [matrixonly]: example = 1, a = 0, b = 1, d = 0.25.             \n Reference: G. M. Wing and J. D. Zahrt, A Primer on Integral 
             Equations of the First Kind, SIAM, Philadelphia, 1991, p. 17.",

             "blur" => "Image deblurring test problem:
             The generated matrix A is an n*n-by-n*n sparse, symmetric, 
             doubly block Toeplitz matrix that models blurring of an n-by-n 
             image by a Gaussian point spread function.
             \n Input options:
             \n 1. [type,] n, band, σ, [matrixonly]: the dimension of the matrix
             is n^2. band is the half-bandwidth, only matrix elements within
             a distance band-1 from the diagonal are nonzero. σ controls the
             width of the Gaussin point spread function. The larger the σ, the 
             wider the function and the more ill posed the problem. 
             If matrixonly = false, the linear system A, b, x will be generated.
             (matrixonly = true by default.)
             \n 2. [type,] n, [matrixonly]: band = 3, σ = 0.7.",
             
             "spikes" => "Test problem with \"spike\" solution:
             Artifically generated discrete ill-posed problem.
             \n Input options:
             \n 1. [type,] n, t_max, [matrixonly]: the dimension of the 
             matrix is n. t_max controls the length of the pulse train.
             If matrixonly = false, the linear system A, b, x will be 
             generated. (matrixonly = true by default.) The solution x
             consists a unit step at t = .5 and a pulse train of spike
             of decreasing magnitude at t = .5, 1.5, 2.5, ...;
             \n 2. [type,] n, [matrixonly]: t_max = 5.    
             ",
             
             "ursell" => "Integral equation with no square integrable 
             solution. Discretization of a first kind Fredholm integral
             equation with kernel K and right-hand side g given by
                     K(s,t) = 1/(s+t+1), g(s) = 1,
             where both integration intervals are [0, 1]. The matrix A
             is a Hankel matrix.
             \n Input options:
             \n 1. [type,] n, [matrixonly]: the dimension of the matrix
             is n. If matrixonly = false, the right-hand side b will also
             be generated. (matrixonly = true by default).
             \n Reference: F. Ursell, Introduction to the theory of linear
             integral equations., Chapter 1 in L. M. Delves & J. Walsh,
             Numerical Solution of Integral Equations, Clarendon Press, 
             1974."
             );

matrixclass =
@compat Dict("symmetric" => ["hilb", "cauchy", "circul", "dingdong",
                             "invhilb", "moler", "pascal", "pei",
                             "clement", "fiedler", "minij", "tridiag",
                             "lehmer", "randcorr", "poisson", "wilkinson",
                              "kms", "wathen", "oscillate", "prolate", 
                             "hankel"],

             "inverse" => ["hilb", "hadamard", "cauchy", "invhilb",
                           "forsythe", "magic", "triw", "moler", "pascal",
                           "kahan", "pei", "vand", "invol", "lotkin",
                           "clement", "fiedler", "minij", "tridiag",
                           "lehmer", "poisson", "kms" ],

             "ill-cond" => ["hilb", "cauchy", "frank", "invhilb",
                            "forsythe", "triw", "moler", "pascal",
                            "kahan","pei", "vand", "invol", "lotkin",
                            "tridiag", "rosser", "randsvd", "kms", 
                            "oscillate", "prolate", "golub"],

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
                          "wathen", "oscillate", "golub"],

             "regprob" => ["deriv2", "shaw", "wing", "foxgood", "heat", 
                           "baart", "phillips", "gravity", "blur", 
                           "spikes", "ursell"],
               );
