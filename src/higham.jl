# Higham Test Matrices

"""
Hilbert matrix
==============
The Hilbert matrix is a very ill conditioned matrix. 
It is symmetric positive definite and totally positive. 

*Input options:*

1. [type,] dim: the dimension of the matrix;

2. [type,] row_dim, col_dim: the row and column dimension.

['inverse', 'ill-cond', 'symmetric', 'pos-def']

*Reference:* M. D. Choi, Tricks or treats with the Hilbert matrix,
Amer. Math. Monthly, 90 (1983), pp. 301-312.
N. J. Higham, Accuracy and Stability of Numerical Algorithms,
Society for Industrial and Applied Mathematics, Philadelphia, PA,
    USA, 2002; sec. 28.1.
"""
function hilb{T}(::Type{T}, m::Int, n::Int)
    # compute the Hilbert matrix
    H = zeros(T, m, n)
    for j = 1:n, i= 1:m
        @inbounds H[i,j] = one(T)/ (i + j - one(T))
    end
    return H
end
hilb{T}(::Type{T}, n::Int) = hilb(T, n, n)
hilb(args...) = hilb(Float64, args...)

"""
Inverse of Hilbert matrix
=========================
*Input options:*

1. [type,] dim: the dimension of the matrix.

['inverse', 'ill-cond', 'symmetric','pos-def']

*Reference:* M. D. Choi, Tricks or treats with the Hilbert matrix,
    Amer. Math. Monthly, 90 (1983), pp. 301-312.
    N. J. Higham, Accuracy and Stability of Numerical Algorithms,
    Society for Industrial and Applied Mathematics, Philadelphia, PA,
    USA, 2002; sec. 28.1.
"""
function invhilb{T}(::Type{T}, n::Int)
    # compute the inverse of Hilbert matrix
    # Type: data type
    # n: the dimension of the matrix
    Inv = zeros(T, n, n)
    for i = 1:n, j = 1:n
        Inv[i,j] = (-1)^(i+j)*(i+j-1)*binomial(n+i-1, n-j)*
        binomial(n+j-1,n-i)*binomial(i+j-2,i-1)^2
    end
    return Inv
end
invhilb(n::Int) = invhilb(Float64, n)

"""
Hadamard matrix:
        The Hadamard matrix is a square matrix whose entries are 
        1 or -1. It was named after Jacques Hadamard. The rows of 
        a Hadamard matrix are orthogonal.
         Input options:
         1. [type,] dim: the dimension of the matrix, n is a power of 2.
         ['inverse', 'orthogonal', 'eigen']
         Reference: S. W. Golomb and L. D. Baumert, The search for
           Hadamard matrices, Amer. Math. Monthly, 70 (1963) pp. 12-17
"""
function hadamard{T}(::Type{T}, n::Int)
    #Compute the hadamard matrix
    if n < 1
        lgn = 0
    else
        lgn = @compat round(Integer, log2(n))
    end
    2^lgn != n && throw(ArgumentError("n must be positive integer and a power of 2."))

    H = reshape(T[1], 1, 1)
    for i = 1:lgn
        H = [[H H]; [H -H]]
    end
    return H
end


"""
Cauchy matrix:
            Given two vectors x and y, the (i,j) entry of the Cauchy 
            matrix is 1/(x[i]-y[j]). 
             Input options:
             1. [type,] vec1, vec2: two vectors;
             2. [type,] vec: a vector;
             3. [type,] dim: the dimension of the matrix.
             ['inverse', 'ill-cond', 'symmetric', 'pos-def']
             Reference:  N. J. Higham, Accuracy and Stability of Numerical Algorithms,
            Society for Industrial and Applied Mathematics, Philadelphia, PA,
                USA, 2002; sec. 28.1
"""
function cauchy{T}(::Type{T}, x::Vector, y::Vector)
    # Compute the cauchy matrix
    m = size(x,1)
    n = size(y,1)
    C = zeros(T, m, n)
    [C[i,j]= one(T)/(x[i] + y[j]) for i = 1:m, j = 1:n]
    return C
end

cauchy{T}(::Type{T}, x::Vector) = cauchy(T, x, x)
cauchy{T}(::Type{T}, k::Int) = cauchy(T, [1:k;])

"""
Circulant matrix:
                A circulant matrix has the property that each row is obtained
                by cyclically permuting the entries of the previous row one
                step forward.
                 Input options:
                 1. [type,] vec, n: a vector and the column dimension;
                 2. [type,] vec: a vector;
                 3. [type,] dim: the dimension of the matrix.
                 ['symmetric', 'pos-def', 'eigen']
                 Reference:  P. J. Davis, Circulant Matrices, John Wiley, 1977.
"""
function circul{T}(::Type{T}, v::Vector, w::Vector)
    # Compute the circul matrix
    # v: the first row of the matrix
    # w: the length of the vector is the dimension of the column
    m = length(v)
    n = length(w)
    C = zeros(T, m, n)
    [C[i,j] = v[1 + mod(j - i, n)] for i = 1:m, j = 1:n]
    return C
end

circul{T}(::Type{T}, v::Vector, n::Int) = circul(T, v, T[1:n;])
circul{T}(::Type{T}, v::Vector) = circul(T, v, v)
circul{T}(::Type{T}, k::Int) = circul(T, [1:k;])

"""
Dingdong matrix:
             The Dingdong matrix is a symmetric Hankel matrix invented
             by DR. F. N. Ris of IBM, Thomas J Watson Research Centre.
             The eigenvalues cluster around π/2 and -π/2.
              Input options:
              1. [type,] n: the dimension of the matrix.
              ['symmetric', 'eigen']
              Reference: J. C. Nash, Compact Numerical Methods for
             Computers: Linear Algebra and Function Minimisation,
             second edition, Adam Hilger, Bristol, 1990 (Appendix 1).
"""
function dingdong{T}(::Type{T}, n::Int)
    # Compute the dingdong matrix
    # type: data type
    # n: the dimension of the matrix
    D = zeros(T, n, n)
    [D[i,j] = one(T)/(2*(n-i-j+1.5)) for i = 1:n, j= 1:n]
    return D
end

"""
Frank matrix:
             The Frank matrix is an upper Hessenberg matrix with 
             determinant 1. The eigenvalues are real, positive and 
             very ill conditioned.
              Input options:
              1. [type,] n, k: n is the dimension of the matrix, k = 0 or 1.
             If k = 1 the matrix reflect about the anti-diagonal;
              2. [type,] n: n is the dimension of the matrix.
              ['ill-cond', 'eigen']
              Reference:  W. L. Frank, Computing eigenvalues of complex matrices
             by determinant evaluation and by methods of Danilewski and Wielandt,
             J. Soc. Indust. Appl. Math., 6 (1958), pp. 378-392 (see pp. 385, 388).
"""
function frank{T}(::Type{T}, n::Int, k::Int)
    # Compute the frank matrix
    # type: data type
    # n: the dimension of the matrix
    # k = 0 or 1: k = 1 reflect about the anti-diagonal
    p = T[n:-1:1;];
    F = triu(ones(T, n, 1)*p') + diagm(p[2:n], -1);
    k == 0 ?  F :
    k == 1 ?  F[p,p]' : error("k = 0 or 1, but get ", k)
end
frank{T}(::Type{T}, n::Int) = frank(T, n, 0)

#
# Jordan block
#
function jordbloc{T}( n::Int, lambda::T)
    # Generate a n-by-n jordan block with eigenvalue
    # lambda
    J = lambda*eye(T, n) + diagm(ones(T, n-1, 1)[:,1] , 1)
    return J
end

"""
Forsythe matrix:
             The Forsythe matrix is a n-by-n perturbed Jordan block.
              Input options:
              1. [type,] n, alpha, lambda: n is the dimension of the matrix.
             alpha and lambda are scalars;
              2. [type,] n: alpha = sqrt(eps(type)) and lambda = 0.
              ['inverse', 'ill-cond', 'eigen']
              This generator is adapted from N. J. Higham's Test Matrix Toolbox.
"""
function forsythe{T}(::Type{T}, n::Int , alpha, lambda)
    # Generate the forsythe matrix
    alpha = convert(T, alpha)
    lambda = convert(T, lambda)
    F = jordbloc(n, lambda);
    F[n,1] = alpha;
    return F
end
forsythe{T}(::Type{T}, n::Int) = forsythe(T, n, sqrt(eps(T)), zero(T))


function oddmagic{T}(::Type{T}, n::Int)
    # compute the magic square of odd orders
    A = zeros(T, n, n)
    i = 1
    j = div(n+1, 2)
    for k = 1:n^2
        is = i
        js = j
        A[i,j] = k
        i = n - rem(n+1-i,n)
        j = rem(j,n) + 1
        if A[i,j] != 0
            i = rem(is,n) + 1
            j = js
        end
    end
    return A
end

"""
Magic square matrix:
             The magic matrix is a matrix with integer entries such that 
             the row elements, column elements, diagonal elements and 
             anti-diagonal elements all add up to the same number.
              Input options:
              1. [type,] dim: the dimension of the matrix.
              ['inverse']
"""
function magic{T}(::Type{T}, n::Int)
    # Compute a magic square of order n
    # Learnt from Cleve Moler, Experiments with MATLAB, 2011
    if mod(n, 2) == 1
        # n is odd
        M = oddmagic(T, n)
    elseif mod(n, 4) == 0
        # n is doubly even
        a = @compat floor(Integer, mod([1:n;], 4)/2)
        B = broadcast(==, a', a)
        M = broadcast(+, T[1:n:n^2;]',T[0:n-1;])
        for i = 1:n, j = 1:n
            B[i,j] == 1 ? M[i,j] = n^2 + one(T) - M[i,j] : M[i,j]
        end
    else
        # n is singly even
        p = div(n,2)
        M = oddmagic(T, p)
        M = [M M+2*p^2; M+3*p^2 M+p^2]
        if n == 2
            return M
        end
        i = [1:p;]
        k = div(n-2, 4)
        j = [[1:k;]; [(n-k+2) : n;]]
        M[[i;i+p],j] = M[[i+p;i],j]
        i = k+1
        j = [1, i]
        M[[i;i+p],j] = M[[i+p;i],j]
    end
    return M
end

"""
Grcar Matrix:
             The Grcar matrix is a Toeplitz matrix with sensitive 
             eigenvalues.
              Input options:
              1. [type,] dim, k: dim is the dimension of the matrix and
             k is the number of superdiagonals;
              2. [type,] dim: the dimension of the matrix.
              ['eigen']
              Reference: J. F. Grcar, Operator coefficient methods
             for linear equations, Report SAND89-8691, Sandia National
             Laboratories, Albuquerque, New Mexico, 1989 (Appendix 2).
"""
function grcar{T}(::Type{T}, n::Int, k::Int = 3)
    # Compute grcar matrix
    G = tril(triu(ones(T, n,n)), k) - diagm(ones(T, n-1), -1)
    return G
end

"""
Triw Matrix:
             Upper triangular matrices discussed by Wilkinson and others.
              Input options:
               1. [type,] row_dim, col_dim, α, k: row_dim and col_dim
             are row and column dimension of the matrix. α is a
             scalar representing the entries on the superdiagonals.
             k is the number of superdiagonals;
              2. [type,] dim: the dimension of the matrix.
              ['inverse', 'ill-cond']
              Reference:  G. H. Golub and J. H. Wilkinson, Ill-conditioned
             eigensystems and the computation of the Jordan canonical form,
             SIAM Review, 18(4), 1976, pp. 578-6
"""
function triw{T}(::Type{T}, m::Int, n::Int, alpha, k::Int)
    alpha = convert(T, alpha)
    A = tril(eye(T, m,n) + alpha * triu(ones(T, m,n), 1), k)
    return A
end
triw{T}(::Type{T}, n::Int) = triw(T, n, n,  -1, n-1)

"""
Moler Matrix:
             The Moler matrix is a symmetric positive definite matrix. 
             It has one small eigenvalue.
              Input options:
              1. [type,] dim, alpha: dim is the dimension of the matrix,
             alpha is a scalar;
              2. [type,] dim: alpha = -1.
              ['inverse', 'ill-cond', 'symmetric', 'pos-def']
              Reference:  J.C. Nash, Compact Numerical Methods for Computers:
             Linear Algebra and Function Minimisation, second edition,
             Adam Hilger, Bristol, 1990 (Appendix 1).
"""
function moler{T}(::Type{T}, n::Int, alpha = -1.)
    alpha = convert(T, alpha)
    M = triw(T, n, n, alpha, n - 1)' * triw(T, n, n, alpha, n -1 )
    return M
end

"""
Pascal Matrix:
             The Pascal matrix’s anti-diagonals form the Pascal’s triangle.
              Input options:
              1. [type,] dim: the dimension of the matrix.
              ['Inverse', 'ill-cond', 'symmetric', 'pos-def', 'eigen']
              Reference: R. Brawer and M. Pirovino, The linear algebra of
             the Pascal matrix, Linear Algebra and Appl., 174 (1992),
             pp. 13-23 (this paper gives a factorization of L = PASCAL(N,1)
             and a formula for the elements of L^k).
             N. J. Higham, Accuracy and Stability of Numerical Algorithms
             Society for Industrial and Applied Mathematics, Philadelphia, PA,
             USA, 2002; sec. 28.4.
"""
function pascal{T}(::Type{T}, n::Int)
    P = zeros(T, n,n)
    for j = 1:n
        for i = 1:n
            try P[i,j] = binomial(i+j-2, j-1) 
            catch 
                P[i,j] = binomial(BigInt(i+j-2), j-1)
            end
        end
    end
    return P
end

"""
Kahan Matrix:
             The Kahan matrix is an upper trapezoidal matrix, i.e., the 
             (i,j) element is equal to 0 if i > j. The useful range of 
             θ is 0 < θ < π. The diagonal is perturbed by 
             pert*eps()*diagm([n:-1:1;]).
              Input options:
              1. [type,] m, n, θ, pert: m, n are the row and column
             dimensions of the matrix. theta and pert are scalars;
              2. [type,] dim, θ, pert: dim is the dimension of the matrix;
              3. [type,] dim: θ = 1.2, pert = 25.
              ['inverse', 'ill-cond']
              Reference: W. Kahan, Numerical linear algebra, Canadian Math.
             Bulletin, 9 (1966), pp. 757-801.
"""
function kahan{T}(::Type{T}, m::Int, n::Int, theta, pert)
    theta = convert(T, theta)
    pert = convert(T, pert)
    s = sin(theta)
    c = cos(theta)
    dim = min(m,n)
    A = zeros(T, m, n)
    for i = 1:m, j = 1:n
        i > dim ? A[i,j] = zero(T) :
        i > j   ? A[i,j] = zero(T) :
        i==j    ? A[i,j] = s^(i-1)+pert*eps(T)*(m-i+1) : A[i,j] = -c*s^(i-1)
    end
    return A
end
kahan{T}(::Type{T}, n::Int, theta, pert) = kahan(T, n, n, theta, pert)
kahan{T}(::Type{T}, n::Int) = kahan(T, n, n, 1.2, 25.)

"""
Pei Matrix:
             The Pei matrix is a symmetric matrix with known inversion.
              Input options:
              1. [type,] dim, alpha: dim is the dimension of the matrix.
             alpha is a scalar;
              2. [type,] dim: the dimension of the matrix.
              ['inverse', 'ill-cond', 'symmetric', 'pos-def']
              Reference: M. L. Pei, A test matrix for inversion procedures,
             Comm. ACM, 5 (1962), p. 508.
"""
function pei{T}(::Type{T}, n::Int, alpha = 1)
    alpha = convert(T, alpha)
    return alpha*eye(T, n, n) + ones(T, n, n)
end

"""
Vandermonde Matrix:
             The inverse and determinat are known explicity.
              Input options:
              1. [type,] vec, dim: vec is a vector, dim is the number of columns.
              2. [type,] vec
              3. [type,] dim
              ['inverse', 'ill-cond']
              Reference: N. J. Higham, Stability analysis of algorithms
             for solving confluent Vandermonde-like systems, SIAM J.
             Matrix Anal. Appl., 11 (1990), pp. 23-41.
"""
function vand{T}(::Type{T}, p::Vector, n::Int)
    # n: number of rows
    # p: a vector
    m = length(p)
    V = Array(T, m, n)
    for j = 1:m
        @inbounds V[j, 1] = 1
    end
    for i = 2:n
        for j = 1:m
            @inbounds V[j,i] = p[j] * V[j,i-1]
        end
    end
    return V
end
vand{T}(::Type{T}, n::Int) = vand(T, [1:n;], n)
vand{T}(::Type{T}, p::Vector) = vand(T, p, length(p))

"""
Involutory Matrix:
             An involutory matrix is a matrix that is its own inverse.
              Input options:
              1. [type,] dim: dim is the dimension of the matrix.
              ['inverse', 'ill-cond', 'eigen']
              Reference: A. S. Householder and J. A. Carpenter, The
             singular values of involutory and of idempotent matrices,
             Numer. Math. 5 (1963), pp. 234-237.
"""
function invol{T}(::Type{T}, n::Int)
    A = hilb(T, n)
    d = -n
    A[:,1] = d*A[:, 1]
    for i = 1:n-1
        d = -(n+i)*(n-i)*d/(i*i)
        A[i+1,:] = d*A[i+1,:]
    end
    return A
end

"""
Chebyshev spectral differentiation matrix:
             If k = 0,the generated matrix is nilpotent and a vector with 
             all one entries is a null vector. If k = 1, the generated 
             matrix is nonsingular and well-conditioned. Its eigenvalues 
             have negative real parts.
              Input options:
              1. [type,] dim, k: dim is the dimension of the matrix and
             k = 0 or 1.
              2. [type,] dim
              ['eigen']
              Reference: L. N. Trefethen and M. R. Trummer, An instability
             phenomenon in spectral methods, SIAM J. Numer. Anal., 24 (1987), pp. 1008-1023.
"""
function chebspec{T}(::Type{T}, n::Int, k::Int = 0)
    # k = 0 or 1
    k == 1 ? n = n + 1 : none
    c = ones(T, n)
    c[1] = 2.
    c[2:n-1] = one(T)
    c[n] = 2.
    x = ones(T, n)

    A = zeros(T, n, n)
    # Compute the chebyshev points
    for i = 1:n
        x[i] = cos(pi * (i - 1) / (n - 1))
    end

    for j = 1:n, i = 1:n
        i != j ? A[i,j] = (-1)^(i+j) * c[i] / (c[j] * (x[i] - x[j])) :
        i == 1 ? A[i,i] = (2 * (n -1)^2 + 1) / 6 :
        i == n ? A[i,i] = - (2 * (n-1)^2 + 1) / 6 :
        A[i,i] = - 0.5 * x[i] / (1 - x[i]^2)
    end
    k == 1 ? A = A[2:n, 2:n] : none
    return A
end

"""
Lotkin Matrix:
             The Lotkin matrix is the Hilbert matrix with its first row 
             altered to all ones. It is unsymmetric, ill-conditioned and 
             has many negative eigenvalues of small magnitude.
              Input options:
              1. [type,] dim: dim is the dimension of the matrix.
              ['inverse', 'ill-cond', 'eigen']
              Reference: M. Lotkin, A set of test matrices, MTAC, 9 (1955), pp. 153-161.
"""
function lotkin{T}(::Type{T}, n::Int)
    A = hilb(T, n)
    A[1,:] = ones(T,n)'
    return A
end

"""
Clement Matrix:
             The Clement matrix is a tridiagonal matrix with zero
             diagonal entries. If k = 1, the matrix is symmetric.
              Input options:
              1. [type,] dim, k: dim is the dimension of the matrix.
             If k = 0, the matrix is Tridiagonal. If k = 1, the matrix
             is SymTridiagonal;
              2. [type,] dim: k = 0.
              ['inverse', 'symmetric', 'eigen']
              Reference: P. A. Clement, A class of triple-diagonal
             matrices for test purposes, SIAM Review, 1 (1959), pp. 50-52.
"""
function clement{T}(::Type{T}, n::Int, k::Int = 0)
    # construct Tridiagonal matrix
    # n is the dimension of the matrix
    # k = 0 or 1
    if n == 1 # handle the 1-d case.
        return zeros(T, 1, 1)
    end
    n = n -1
    x = T[n:-1:1;]
    z = T[1:n;]
    if k == 0
        A = Tridiagonal(x,zeros(T,n+1),z)
    else
        y = sqrt(x.*z)
        A = SymTridiagonal(zeros(T,n+1), y)
    end
    return A
end

"""
Fiedler Matrix:
             The Fiedler matrix is symmetric matrix with a dominant 
             positive eigenvalue and all the other eigenvalues are negative.
              Input options:
              1. [type,] vec: vec is a vector;
              2. [type,] dim: dim is the dimension of the matrix.
              ['inverse', 'symmetric', 'eigen']
              Reference: G. Szego, Solution to problem 3705, Amer. Math.
              Monthly, 43 (1936), pp. 246-259.
             J. Todd, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
             Birkhauser, Basel, and Academic Press, New York, 1977, p. 159.
"""
function fiedler{T}(::Type{T}, v::Vector)
    n = length(v)
    v = v[:].'
    A = ones(T, n) * v
    A = abs(A - A.') # nonconjugate transpose
end
fiedler{T}(::Type{T}, n::Int) = fiedler(T, [1:n;])

"""
MIN[I,J] Matrix:
             A matrix with (i,j) entry min(i,j). It is a symmetric positive
             definite matrix. The eigenvalues and eigenvectors are known 
             explicitly. Its inverse is tridiagonal.
              Input options:
              1. [type,] dim: dim is the dimension of the matrix.
              ['inverse', 'symmetric', 'pos-def', 'eigen']
              Reference: J. Fortiana and C. M. Cuadras, A family of matrices,
             the discretized Brownian bridge, and distance-based regression,
             Linear Algebra Appl., 264 (1997), 173-188.  (For the eigensystem of A.)
"""
function minij{T}(::Type{T}, n::Int)
    A = zeros(T, n, n)
    [A[i,j] = min(i,j) for i = 1:n, j = 1:n]
    return A
end

"""
Binomial Matrix:
             The matrix is a multiple of an involutory matrix.
              Input options:
              1. [type,] dim: dim is the dimension of the matrix.
              ['']
"""
function binomialm{T}(::Type{T}, n::Int)
    # Mulitiple of involutory matrix
    L = Array(T, n, n)
    D = Diagonal((-2).^[0:n-1;])
    [L[i,j] = binomial(i-1, j-1) for i = 1:n, j = 1:n]
    U = L[n:-1:1, n:-1:1]
    return L*D*U
end

"""
Tridiagonal Matrix:
              Input options:
              1. [type,] v1, v2, v3: v1 and v3 are vectors of subdiagonal
             and superdiagonal elements, respectively, and v2 is a vector
             of diagonal elements;
              2. [type,] dim, x, y, z: dim is the dimension of the matrix,
             x, y, z are scalars. x and z are the subdiagonal and 
             superdiagonal elements, respectively, and y is the diagonal 
             elements;
              3. [type,] dim: x = -1, y = 2, z = -1. This matrix is also
             known as the second difference matrix.
              ['inverse', 'ill-cond', 'pos-def', 'eigen']
              Reference: J. Todd, Basic Numerical Mathematics, Vol. 2:
             Numerical Algebra, Birkhauser, Basel, and Academic Press,
             New York, 1977, p. 155.
"""
function tridiag{T}(::Type{T}, x::Vector, y::Vector, z::Vector)
    x = map((i)-> convert(T, i), x)
    y = map((i)-> convert(T, i), y)
    z = map((i)-> convert(T, i), z)
    return Tridiagonal(x,y,z)
end
# Toeplitz tridiagonal matrix
tridiag{T}(::Type{T}, n::Int, x::Int, y::Int, z::Int) =
n == 1 ? y*ones(T,1,1) :
         tridiag(T, x*ones(T, n-1), y*ones(T, n), z*ones(T, n-1))
tridiag{T}(::Type{T}, n::Int) = tridiag(T, n, -1, 2, -1)

"""
Lehmer Matrix:
             The Lehmer matrix is a symmetric positive definite matrix. 
             It is totally nonnegative. The inverse is tridiagonal and 
             explicitly known
              Input options:
              1. [type,] dim: the dimension of the matrix.
              ['inverse', 'symmetric', 'pos-def']
              Reference: M. Newman and J. Todd, The evaluation of
             matrix inversion programs, J. Soc. Indust. Appl. Math.,
             6 (1958), pp. 466-476.
             Solutions to problem E710 (proposed by D.H. Lehmer): The inverse
             of a matrix, Amer. Math. Monthly, 53 (1946), pp. 534-535.
"""
function lehmer{T}(::Type{T}, n::Int)
    A = Array(T, n, n)
    [A[i,j] = min(i,j) / max(i,j) for i = 1:n, j = 1:n]
    return A
end

"""
Parter Matrix:
             The Parter matrix is a Toeplitz and Cauchy matrix 
             with singular values near π.
              Input options:
              1. [type,] dim: the dimension of the matrix.
              ['eigen']
              Reference: The MathWorks Newsletter, Volume 1, Issue 1,
             March 1986, page 2. S. V. Parter, On the distribution of the
             singular values of Toeplitz matrices, Linear Algebra and
             Appl., 80 (1986), pp. 115-130.
"""
function parter{T}(::Type{T}, n::Int)
    A = Array(T, n, n)
    [A[i,j] = one(T) / (i - j + 0.5) for i = 1:n, j = 1:n]
    return A
end

"""
Chow Matrix:
             The Chow matrix is a singular Toeplitz lower Hessenberg matrix.
              Input options:
              1. [type,] dim, alpha, delta: dim is dimension of the matrix.
             alpha, delta are scalars such that A[i,i] = alpha + delta and
             A[i,j] = alpha^(i + 1 -j) for j + 1 <= i;
              2. dim: alpha = 1, delta = 0.
              ['eigen']
              Reference: T. S. Chow, A class of Hessenberg matrices with known
             eigenvalues and inverses, SIAM Review, 11 (1969), pp. 391-395.
"""
function chow{T}(::Type{T}, n::Int, alpha, delta)
    A = zeros(T, n, n)
    alpha = convert(T, alpha)
    delta = convert(T, delta)
    for i = 1:n, j = 1:n
        if i == j - 1
            A[i,j] = one(T)
        elseif i == j
            A[i,j] = alpha + delta
        elseif j + 1 <= i
            A[i,j] = alpha^(i + 1 - j)
        end
    end
    return A
end
chow{T}(::Type{T}, n::Int) = chow(T, n, 1, 0)

#
# newsign: newsign(0) = 1
#
function newsign(x)
    x == 0 ? y = 1 : y = sign(x)
    return y
end

"""
Random Correlation Matrix:
             A random correlation matrix is a symmetric positive 
             semidefinite matrix with 1s on the diagonal.
              Input options:
              1. [type,] dim: the dimension of the matrix.
              ['symmetric', 'pos-semidef', 'random']
"""
function randcorr{T}(::Type{T}, n::Int)
    x = rand(T,n) # x is the vector of random eigenvalues from a uniform distribution.
    x = n * x / sum(x) # x has nonnegtive elements.
    A = diagm(x)
    F = qrfact(randn(n,n));
    Q = F[:Q]*diagm(sign(diag(F[:R]))) # form a random orthogonal matrix.
    copy!(A, Q*A*Q')

    a = diag(A)
    l = find(a .< 1)
    g = find(a .> 1)

    # Apply Given rotation to set A[i,i] = 1
    while length(l) > 0 && length(g) > 0
        k =  @compat ceil(Integer, rand()*length(l))
        h =  @compat ceil(Integer, rand()*length(g))
        i = l[k]
        j = g[h]
        if i > j
            i,j = j,i
        end
        alpha = sqrt(A[i,j]^2 - (a[i] - 1)*(a[j] - 1))
        # take sign to avoid cancellation.
        t = (A[i,j] + newsign(A[i,j]) * alpha) / (a[j] - 1)
        c = 1/ sqrt(1 + t^2)
        s = t*c

        A[:, [i,j]] = A[:, [i,j]] * [c s; -s c]
        A[[i,j], :] = [c -s; s c] * A[[i,j], :]

        A[i,i] = 1
        a = diag(A)
        l = find(a.<1)
        g = find(a.>1)
    end
    [A[i,i] = 1 for i = 1:n]
    return (A + A')/2

end

"""
Poisson Matrix:
             A block tridiagonal matrix from Poisson’s equation. 
             This matrix is sparse, symmetric positive definite and 
             has known eigenvalues.
              Input options:
              1. [type,] n: the dimension of the matirx is n^2.
              ['inverse', 'symmetric', 'pos-def', 'eigen', 'sparse']
              Reference: G. H. Golub and C. F. Van Loan, Matrix Computations,
             second edition, Johns Hopkins University Press, Baltimore,
             Maryland, 1989 (Section 4.5.4).
"""
function poisson{T}(::Type{T}, n::Int)
    S = full(tridiag(T, n))
    A = speye(T, n)
    return kron(A,S) + kron(S,A)
end

#
# Toeplitz Matrix
#
function toeplitz{T}(::Type{T}, vc::Vector, vr::Vector)
    n = length(vc)
    length(vr) == n || throw(DimensionMismatch(""))
    vc[1] == vr[1] || error("The first element of the vectors must be the same.")    
    A = Array(T, n, n)
    [i>=j ? A[i,j] = vc[i-j+1]: A[i,j] = vr[j-i+1] for i=1:n, j=1:n]
    A
end
toeplitz{T}(::Type{T}, v::Vector) = toeplitz(T, v, v)
toeplitz{T}(::Type{T}, n::Int) = toeplitz(T, [1:n;])


#
# Hankel matrix
#
function hankel{T}(::Type{T}, vc::Vector, vr::Vector)
    p = [vc; vr[2:end]]
    m = length(vc)
    n = length(vr)
    H = Array(T, m, n)
    [H[i,j] = p[i+j-1] for i=1:m, j=1:n]
    H
end
hankel{T}(::Type{T}, v::Vector) = hankel(T, v, v)
hankel{T}(::Type{T}, n::Int) = hankel(T, [1:n;])

#
# Prolate Matrix
#
function prolate{T}(::Type{T}, n::Int, w::Real)
    v = Array(T, n)
    v[1] = 2*w
    [v[i] = sin(2*pi*w*i)/pi*i for i = 2:n]
    return toeplitz(T, v)
end
prolate{T}(::Type{T}, n::Int) = prolate(T, n, 0.25)

"""
Neumann Matrix:
             A singular matrix from the discrete Neumann problem. 
             This matrix is sparse and the null space is formed by 
             a vector of ones
              Input options:
              1. [type,] n: the dimension of the matrix is n^2.
              ['eigen', 'sparse']
              Reference: R. J. Plemmons, Regular splittings and the
             discrete Neumann problem, Numer. Math., 25 (1976), pp. 153-161.
"""
function neumann{T}(::Type{T}, n::Int)
    if n == 1
        return 4 * ones(T, 1,1) #handle 1-d case.
    end
    S = full(tridiag(T, n))
    S[1,2] = -2
    S[n, n-1] = -2
    A = speye(T, n)
    return kron(S,A) + kron(A,S)
end

#
# Sylvester's orthogonal matrix
# See Rosser matrix Reference 2.
#
# for a = d = 2, b = c = 1, P_block' * P_block = 10 * Identity
#
P_block{T}(::Type{T}, a, b, c, d) = reshape(T[a, b, c, d, b, -a, -d, c, c, d, -a, -b, d, -c, b, -a], 4,4)

"""
Rosser Matrix:
             The Rosser matrix’s eigenvalues are very close together 
             so it is a challenging matrix for many eigenvalue algorithms.
              Input options:
              1. [type,] dim, a, b: dim is the dimension of the matrix.
             dim must be a power of 2.
             a and b are scalars. For dim = 8, a = 2 and b = 1, the generated
             matrix is the test matrix used by Rosser.
              2. [type,] dim: a = rand(1:5), b = rand(1:5).
              ['eigen', 'ill-cond', 'random']
              Reference: J. B. Rosser, C. Lanczos, M. R. Hestenes, W. Karush,
             Separation of close eigenvalues of a real symmetric matrix,
             Journal of Research of the National Bureau of Standards, v(47)
             (1951)
"""
function rosser{T}(::Type{T}, n::Int, a, b)
    if n < 1
        lgn = 0
    else
        lgn = @compat round(Integer, log2(n))
    end
    2^lgn != n && throw(ArgumentError("n must be positive integer and a power of 2."))

    if n == 1 # handle 1-d case
        return 611 * ones(T, 1, 1)
    end

    if n == 2
        #eigenvalues are 500, 510
        B = T[101 1; 1 101]
        P = T[2 1;1 -2]
        A = P'*B*P
    elseif n == 4
        # eigenvalues are 0.1, 1019.9, 1020, 1020 for a = 2 and b = 1
        B = zeros(T, n, n)
        B[1,1], B[1,4], B[4,1], B[4,4] = 101, 1, 1, 101;
        B[2,2], B[2,3], B[3,2], B[3,3] = 1, 10, 10, 101;
        P = P_block(T, a, b, b, a)
        A = P' * B * P
    elseif n == 8
        # eigenvalues are 1020, 1020, 1000, 1000, 0.098, 0, -1020
        B = zeros(T, n, n)
        B[1,1], B[6,1], B[2,2], B[8,2] = 102, 1, 101, 1;
        B[3,3], B[7,3] = 98, 14;
        B[4,4], B[5,4], B[4,5], B[5,5] = 1, 10, 10, 101;
        B[1,6], B[6,6], B[3,7],B[7,7], B[2,8], B[8,8] = 1, -102, 14, 2, 1, 101;
        P = [P_block(T, a, b, b, a)' zeros(T, 4,4); zeros(T, 4,4) P_block(T, b, -b, -a, a)]
        A = P' * B * P
    else
        lgn = lgn - 2
        halfn = @compat round(Integer, n/2)
        # using Sylvester's method
        P = P_block(T, a, b, b, a)
        m = 4
        for i in 1:lgn
            P = [P zeros(T, m, m); zeros(T, m, m) P]
            m = m * 2
        end
        # mix 4 2-by-2 matrices (with close eigenvalues) into a large nxn matrix.
        B_list = T[102, 1, 1, - 102, 101, 1, 1, 101, 1, 10, 10, 101, 98, 14, 14, 2]
        B = zeros(T, n, n)
        j, k = 1, 5
        for i in 1:(halfn + 1)
            indexend = halfn -1 + i
            list_start = k
            list_end = k + 3

            if list_start > 16 || list_end > 16
                k = 1
                list_start = 1
                list_end = 4
            end
            B[j,j], B[j,indexend], B[indexend, j], B[indexend, indexend] = B_list[list_start:list_end]
            j = j + 1
            k = k + 4
        end
        A = P' * B * P
    end

    return A
end
rosser{T}(::Type{T}, n::Int) = rosser(T, n, rand(1:5), rand(1:5))


"""
Matrix with Application in Sampling Theory:
             A nonsymmetric matrix with eigenvalues 0, 1, 2, ... n-1.
              Input options:
              1. [type,] vec: vec is a vector with no repeated elements;
              2. [type,] n: the dimension of the matrix is n.
              ['eigen']
              Reference: L. Bondesson and I. Traat, A nonsymmetric matrix
             with integer eigenvalues, linear and multilinear algebra, 55(3)
             (2007), pp. 239-247
"""
function sampling{T}(::Type{T}, x::Vector)
    n = length(x)
    A = zeros(T, n, n)
    for j = 1:n, i = 1:n
        if i != j
            A[i,j] = x[i] / (x[i] - x[j])
        end
    end
    d = sum(A, 2)
    A = A + diagm(d[:])
    return A
end
#
# special probability case
# see:
#   L. Bondesson and I. Traat, A Nonsymmetric Matrix with Integer
#   Eigenvalues, Linear and Multilinear Algebra, 55(3)(2007), pp. 239-247.
#
function sampling{T}(::Type{T}, n::Int)
    p = T[1:n;] / n
    return sampling(T, p)
end

"""
Wilkinson Matrix:
             The Wilkinson matrix is a symmetric tridiagonal matrix with 
             pairs of nearly equal eigenvalues. The most frequently used                 case is matrixdepot(\"wilkinson\", 21).
              Input options:
              1. [type,] dim: the dimension of the matrix.
              ['symmetric', 'eigen']
              Reference: J. H. Wilkinson, Error analysis of direct methods
             of matrix inversion, J. Assoc. Comput. Mach., 8 (1961),  pp. 281-330.
"""
function wilkinson{T}(::Type{T}, n::Int)
    if n == 1 # handle 1-d case
        return ones(T, 1, 1)
    end
    m = (n-1)/2
    A = Tridiagonal(ones(T,n-1), abs(T[-m:m;]), ones(T, n-1))
    return A
end

"""
Random Matrix with Element -1, 0, 1:
              Input options:
              1. [type,] m, n, k: m and n are row and column dimensions,
             k = 1: entries are 0 or 1.
             k = 2: entries are -1 or 1.
             k = 3: entries are -1, 0 or 1.
              2. [type,] n, k: m = n;
              3. [type,] n: k = 1.
              ['random']
"""
function rando{T}(::Type{T}, m::Int, n::Int, k::Int)
    A = Array(T, m, n)
    if k == 1
        copy!(A, floor(rand(m,n) + .5))
    elseif k == 2
        copy!(A, 2 * floor(rand(m,n) + .5) - one(T))
    elseif k == 3
        copy!(A, round(3 * rand(m,n) - 1.5))
    else
        throw(ArgumentError("invalid k value."))
    end
    return A
end
rando{T}(::Type{T}, n::Int, k::Int) = rando(T, n, n, k)
rando{T}(::Type{T}, n::Int) = rando(T, n, n, 1)

#
# Pre-multiply by random orthogonal matrix
#
function qmult!{T}(A::Matrix{T})
    n, m = size(A)

    d = zeros(T, n)
    for k = n-1:-1:1

        # generate random Householder transformation
        x = randn(n-k+1)
        s = norm(x)
        sgn = sign(x[1]) + (x[1]==0)
        s = sgn * s
        d[k] = -sgn
        x[1] = x[1] + s
        beta = s * x[1]

        # apply the transformation to A
        y = x'*A[k:n, :];
        A[k:n, :] = A[k:n, :] - x * (y /beta)
    end

    # tidy up signs
    for i=1:n-1
        A[i, :] = d[i] * A[i, :]
    end
    A[n, :] = A[n, :] * sign(randn())
    return A
end

"""
Random Matrix with Pre-assigned Singular Values:
              Input options:
              1. [type,] m, n, kappa, mode: m, n are the dimensions of the matrix.
             kappa is the condition number of the matrix.
             mode = 1: one large singular value.
             mode = 2: one small singular value.
             mode = 3: geometrically distributed singular values.
             mode = 4: arithmetrically distributed singular values.
             mode = 5: random singular values with  unif. dist. logarithm;
              2. [type,] n, kappa, mode: m = n;
              3. [type,] n, kappa: mode = 3;
              4. [type,] n: kappa = sqrt(1/eps()), mode = 3.
              ['ill-cond', 'random']
              Reference: N. J. Higham, Accuracy and Stability of Numerical
             Algorithms, Society for Industrial and Applied Mathematics,
             Philadelphia, PA, USA, 2002; sec. 28.3.
"""
function randsvd{T}(::Type{T}, m::Int, n::Int, kappa, mode::Int)
    kappa >= 1 || throw(ArgumentError("Condition number must be at least 1."))
    kappa = convert(T, kappa)

    p = min(m,n)
    if p == 1 # handle 1-d case
        return ones(T, 1, 1)*kappa
    end

    if mode == 3
        factor = kappa^(-1/(p-1))
        sigma = factor.^[0:p-1;]
    elseif mode == 4
        sigma = ones(T, p) - T[0:p-1;]/(p-1)*(1 - 1/kappa)
    elseif mode == 5
        sigma = exp(-rand(p) * log(kappa))
    elseif mode == 2
        sigma = ones(T, p)
        sigma[p] = one(T)/kappa
    elseif mode == 1
        sigma = ones(p)./kappa
        sigma[1] = one(T)
    else
        throw(ArgumentError("invalid mode value."))
    end

    A = zeros(T, m, n)
    A[1:p, 1:p] = diagm(sigma)
    A = qmult!(A')
    A = qmult!(A')

    return A
end
randsvd{T}(::Type{T}, n::Int, kappa, mode) = randsvd(T, n, n, kappa, mode)
randsvd{T}(::Type{T}, n::Int, kappa) = randsvd(T, n, kappa, 3)
randsvd{T}(::Type{T}, n::Int) = randsvd(T, n, sqrt(1/eps(T)))

"""
Random Orthogonal Upper Hessenberg Matrix:
             The matrix is constructed via a product of Givens rotations.
              Input options:
              1. [type,] n : n is the dimension of the matrix.
              ['random']
              Reference:  W. B. Gragg, The QR algorithm for unitary
             Hessenberg matrices, J. Comp. Appl. Math., 16 (1986), pp. 1-8.
"""
function rohess{T}(::Type{T}, n::Int)
    x = rand(n-1)*2*pi
    H = eye(T, n)
    H[n,n] = sign(randn())
    for i = n:-1:2
        theta = x[i-1]
        c = convert(T, cos(theta))
        s = convert(T, sin(theta))
        H[[i-1; i], :] = [c*H[i-1, :][:]' + s*H[i, :][:]'; -s*H[i-1, :][:]' + c*H[i, :][:]']
    end
    return H
end


"""
Kac-Murdock-Szego Toeplitz Matrix:
              Input:
              1. [type,] n, rho: n is the dimension of the matrix, rho is a
             scalar such that A[i,j] = rho^(abs(i-j)).
              2. [type,] n: rho = 0.5
              ['inverse', 'ill-cond', 'symmetric', 'pos-def']
              Reference: W. F. Trench, Numerical solution of the eigenvalue
             problem for Hermitian Toeplitz matrices, SIAM J. Matrix Analysis
             and Appl., 10 (1989), pp. 135-146 (and see the references therein).
"""
function kms{T}(::Type{T}, n::Int, rho::Number)
    typeof(rho) <: Complex ? A = Array(typeof(rho), n, n): A = Array(T, n, n)
    [A[i,j] = rho^(abs(i-j)) for i = 1:n, j = 1:n]
    if typeof(rho) <: Complex
        A = conj(tril(A, -1)) + triu(A)
    end
    return A
end
kms{T}(::Type{T}, n::Int) = kms(T, n, convert(T, 0.5))

"""
Wathen Matrix:
             Wathen Matrix is a sparse, symmetric positive, random matrix 
             arose from the finite element method. The generated matrix is 
             the consistent mass matrix for a regular nx-by-ny grid of 
             8-nodes.
              Input options:
              1. [type,] nx, ny: the dimension of the matrix is equal to
             3 * nx * ny + 2 * nx * ny + 1;
              2. [type,] n: nx = ny = n.
              ['symmetric', 'pos-def', 'eigen', 'random', 'sparse']
              Reference: A. J. Wathen, Realistic eigenvalue bounds for
             the Galerkin mass matrix, IMA J. Numer. Anal., 7 (1987),
             pp. 449-457.
"""
function wathen{T}(::Type{T}, nx::Int, ny::Int)
    e1 = T[6 -6 2 -8;-6 32 -6 20;2 -6 6 -6;-8 20 -6 32]
    e2 = T[3 -8 2 -6;-8 16 -8 20;2 -8 3 -8;-6 20 -8 16]
    e3 = [e1 e2; e2' e1]/45
    n = 3 * nx * ny + 2 * nx + 2 * ny + 1
    ntriplets = nx * ny * 64
    Irow = zeros(Int, ntriplets)
    Jrow = zeros(Int, ntriplets)
    Xrow = zeros(T, ntriplets)
    ntriplets = 0
    rho = 100 * rand(nx, ny)
    node = zeros(T, 8)

    for j = 1:ny
        for i = 1:nx

            node[1] = 3 * j * nx + 2 * i + 2 * j + 1
            node[2] = node[1] - 1
            node[3] = node[2] - 1
            node[4] = (3 * j - 1) * nx + 2 * j + i - 1
            node[5] = (3 * j - 3) * nx + 2 * j + 2 * i - 3
            node[6] = node[5] + 1
            node[7] = node[5] + 2
            node[8] = node[4] + 1

            em = convert(T, rho[i,j]) * e3

            for krow = 1:8
                for kcol = 1:8
                    ntriplets +=  1
                    Irow[ntriplets] = node[krow]
                    Jrow[ntriplets] = node[kcol]
                    Xrow[ntriplets] = em[krow, kcol]
                end
            end

        end
    end
    return sparse(Irow, Jrow, Xrow, n, n)
end
wathen{T}(::Type{T}, n::Int) = wathen(T, n, n)

"""
Golub Matrix:
             Golub matrix is the product of two random unit lower and upper
             triangular matrices respectively. LU factorization without pivoting
             fails to reveal that such matrices are badly conditioned.
              Input options:
              1. [type,] n: the dimension of the matrix is n.
              Reference: D. Viswanath and N. Trefethen. Condition Numbers of 
             Random Triangular Matrices, SIAM J. Matrix Anal. Appl. 19, 564-581,
             1998.
"""
function golub{T}(::Type{T}, n::Int)
    s = 10
    L = Array(T, n, n)
    U = Array(T, n, n)
    if T <: Integer
        [L[i,j] = round_matlab(T, s*randn()) for j = 1:n, i = 1:n]
        [U[i,j] = round_matlab(T, s*randn()) for j = 1:n, i = 1:n]
    else
        [L[i,j] = s*randn() for j = 1:n, i = 1:n]
        [U[i,j] = s*randn() for j = 1:n, i = 1:n]   
    end
    L = tril(L, -1) + eye(L)
    U = triu(U, 1) + eye(U)
    return L*U
end

"""
Companion matrix:
             The companion matrix to a monic polynomial 
                 a(x) = a_0 + a_1x + ... + a_{n-1}x^{n-1} + x^n
             is the n-by-n matrix with ones on the subdiagonal and 
             the last column given by the coefficients of a(x).     
              Input options:
              1. [type,] vec: vec is a vector of coefficients;
              2. [type,] n: vec = [1:n;], the dimension of the matrix 
             is n.
"""
function companion{T}(::Type{T}, v::AbstractVector)
    n = length(v)
    A = zeros(T, n, n)
    A[:, end] = v
    for i = 1:n-1
        A[i+1, i] = one(T)
    end
    A
end
companion{T}(::Type{T}, n::Int) = companion(T, [1:n;])
