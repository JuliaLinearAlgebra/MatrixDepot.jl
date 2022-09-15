# Higham Test Matrices

"""
Hilbert Matrix
==============
The Hilbert matrix has `(i,j)` element `1/(i+j-1)`. It is
notorious for being ill conditioned. It is symmetric
positive definite and totally positive.

*Input options:*

+ [type,] dim: the dimension of the matrix.

+ [type,] row_dim, col_dim: the row and column dimensions.

*Groups:* ["inverse", "illcond", "symmetric", "posdef"]

*References:*

**M. D. Choi**, Tricks or treats with the Hilbert matrix,
Amer. Math. Monthly, 90 (1983), pp. 301-312.

**N. J. Higham**, Accuracy and Stability of Numerical Algorithms,
second edition, Society for Industrial and Applied Mathematics,
Philadelphia, PA, USA, 2002; sec. 28.1.
"""
function hilb(::Type{T}, m::Integer, n::Integer) where {T}
    # compute the Hilbert matrix
    H = zeros(T, m, n)
    for j = 1:n, i = 1:m
        @inbounds H[i, j] = one(T) / (i + j - one(T))
    end
    return H
end
hilb(::Type{T}, n::Integer) where {T} = hilb(T, n, n)
hilb(args...) = hilb(Float64, args...)
hilb(::Type, args...) = throw(MethodError(hilb, Tuple(args)))

"""
Inverse of the Hilbert Matrix
=============================
*Input options:*

+ [type,] dim: the dimension of the matrix.

*Groups:* ["inverse", "illcond", "symmetric","posdef"]

*References:*

**M. D. Choi**, Tricks or treats with the Hilbert matrix,
    Amer. Math. Monthly, 90 (1983), pp. 301-312.

**N. J. Higham**, Accuracy and Stability of Numerical Algorithms, second
    edition, Society for Industrial and Applied Mathematics, Philadelphia, PA,
    USA, 2002; sec. 28.1.
"""
function invhilb(::Type{T}, n::Integer) where {T}
    # compute the inverse of Hilbert matrix
    # Type: data type
    # n: the dimension of the matrix
    Inv = zeros(T, n, n)
    for i = 1:n, j = 1:n
        Inv[i, j] = (-1)^(i + j) * (i + j - 1) * binomial(n + i - 1, n - j) *
                    binomial(n + j - 1, n - i) * binomial(i + j - 2, i - 1)^2
    end
    return Inv
end
invhilb(n::Integer) = invhilb(Float64, n)

"""
Hadamard Matrix
===============
The Hadamard matrix is a square matrix whose entries are
1 or -1. It was named after Jacques Hadamard. The rows of
a Hadamard matrix are orthogonal.

*Input options:*

+ [type,] dim: the dimension of the matrix, `dim` is a power of 2.

*Groups:* ["inverse", "orthogonal", "eigen"]

*References:*

**S. W. Golomb and L. D. Baumert**, The search for
Hadamard matrices, Amer. Math. Monthly, 70 (1963) pp. 12-17
"""
function hadamard(::Type{T}, n::Integer) where {T}
    #Compute the hadamard matrix
    if n < 1
        lgn = 0
    else
        lgn = round(Integer, log2(n))
    end
    2^lgn != n && throw(ArgumentError("n must be positive integer and a power of 2."))

    H = reshape(T[1], 1, 1)
    for i = 1:lgn
        H = [[H H]; [H -H]]
    end
    return H
end
hadamard(n::Integer) = hadamard(Float64, n)

"""
Cauchy Matrix
=============
Given two vectors `x` and `y`, the `(i,j)` entry of the Cauchy matrix is
`1/(x[i]+y[j])`.

*Input options*:

 + [type,] x, y: two vectors.

 + [type,] x: a vector. `y` defaults to `x`.

 + [type,] dim: the dimension of the matrix. `x` and `y` default to
    `[1:dim;]`.

*Groups:* ["inverse", "illcond", "symmetric", "posdef"]

*References:*

**N. J. Higham**, Accuracy and Stability of Numerical Algorithms,
second edition, Society for Industrial and Applied Mathematics, Philadelphia, PA, USA,
2002; sec. 28.1
"""
function cauchy(::Type{T}, x::Vector, y::Vector) where {T}
    # Compute the cauchy matrix
    m = size(x, 1)
    n = size(y, 1)
    C = zeros(T, m, n)
    [C[i, j] = one(T) / (x[i] + y[j]) for i = 1:m, j = 1:n]
    return C
end

cauchy(::Type{T}, x::Vector) where {T} = cauchy(T, x, x)
cauchy(::Type{T}, k::Integer) where {T} = cauchy(T, [1:k;])
cauchy(arg...) = cauchy(Float64, arg...)
cauchy(::Type, args...) = throw(MethodError(cauchy, Tuple(args)))

"""
Circulant Matrix
================
A circulant matrix has the property that each row is obtained
by cyclically permuting the entries of the previous row one
step forward.

*Input options:*

+ [type,] vec, col_dim: a vector and the column dimension.

+ [type,] vec: a vector.

+ [type,] dim: the dimension of the matrix.

*Groups:* ["symmetric", "posdef", "eigen"]

*References:*

**P. J. Davis**, Circulant Matrices, John Wiley, 1977.
"""
function circul(::Type{T}, v::Vector, w::Vector) where {T}
    # Compute the circul matrix
    # v: the first row of the matrix
    # w: the length of the vector is the dimension of the column
    m = length(v)
    n = length(w)
    C = zeros(T, m, n)
    [C[i, j] = v[1+mod(j - i, n)] for i = 1:m, j = 1:n]
    return C
end
circul(::Type{T}, v::Vector, n::Integer) where {T} = circul(T, v, T[1:n;])
circul(::Type{T}, v::Vector) where {T} = circul(T, v, v)
circul(::Type{T}, k::Integer) where {T} = circul(T, [1:k;])
circul(arg...) = circul(Float64, arg...)

"""
Dingdong Matrix
===============
The Dingdong matrix is a symmetric Hankel matrix invented
by DR. F. N. Ris of IBM, Thomas J Watson Research Centre.
The eigenvalues cluster around `π/2` and `-π/2`.

*Input options:*

+ [type,] dim: the dimension of the matrix.

*Groups:* ["symmetric", "eigen"]

*References:*

**J. C. Nash**, Compact Numerical Methods for
Computers: Linear Algebra and Function Minimisation,
second edition, Adam Hilger, Bristol, 1990 (Appendix 1).
"""
function dingdong(::Type{T}, n::Integer) where {T}
    # Compute the dingdong matrix
    # type: data type
    # n: the dimension of the matrix
    D = zeros(T, n, n)
    [D[i, j] = one(T) / (2 * (n - i - j + 1.5)) for i = 1:n, j = 1:n]
    return D
end
dingdong(args...) = dingdong(Float64, args...)
dingdong(::Type, args...) = throw(MethodError(dingdong, Tuple(args)))

"""
Frank Matrix
============
The Frank matrix is an upper Hessenberg matrix with
determinant 1. The eigenvalues are real, positive and
very ill conditioned.

*Input options:*

+ [type,] dim, k: `dim` is the dimension of the matrix, `k = 0 or 1`.
    If `k = 1` the matrix reflect about the anti-diagonal.

+ [type,] dim: the dimension of the matrix.

*Groups:* ["illcond", "eigen"]

*References:*

**W. L. Frank**, Computing eigenvalues of complex matrices
    by determinant evaluation and by methods of Danilewski and Wielandt,
    J. Soc. Indust. Appl. Math., 6 (1958), pp. 378-392 (see pp. 385, 388).
"""
function frank(::Type{T}, n::Integer, k::Integer) where {T}
    # Compute the frank matrix
    # type: data type
    # n: the dimension of the matrix
    # k = 0 or 1: k = 1 reflect about the anti-diagonal
    p = T[n:-1:1;]
    F = triu(ones(T, n, 1) * p') + diagm(-1 => p[2:n])
    k == 0 ? F :
    k == 1 ? F[n:-1:1, n:-1:1]' : error("k = 0 or 1, but get ", k)
end
frank(::Type{T}, n::Integer) where {T} = frank(T, n, 0)
frank(args...) = frank(Float64, args...)
frank(::Type, args...) = throw(MethodError(frank, Tuple(args)))

#
# Jordan block
#
function jordbloc(n::Integer, lambda::T) where {T}
    # Generate a n-by-n jordan block with eigenvalue
    # lambda
    J = lambda * Matrix{T}(I, n, n) + diagm(1 => ones(T, n - 1, 1)[:, 1])
    return J
end

"""
Forsythe Matrix
===============
The Forsythe matrix is a n-by-n perturbed Jordan block.
This generator is adapted from N. J. Higham's Test Matrix Toolbox.

*Input options:*

+ [type,] dim, alpha, lambda: `dim` is the dimension of the matrix.
    `alpha` and `lambda` are scalars.

+ [type,] dim: `alpha = sqrt(eps(type))` and `lambda = 0`.


*Groups:* ["inverse", "illcond", "eigen"]
"""
function forsythe(::Type{T}, n::Integer, alpha, lambda) where {T}
    # Generate the forsythe matrix
    alpha = convert(T, alpha)
    lambda = convert(T, lambda)
    F = jordbloc(n, lambda)
    F[n, 1] = alpha
    return F
end
forsythe(::Type{T}, n::Integer) where {T} = forsythe(T, n, sqrt(eps(T)), zero(T))
forsythe(args...) = forsythe(Float64, args...)
forsythe(::Type, args...) = throw(MethodError(forsythe, Tuple(args)))

function oddmagic(::Type{T}, n::Integer) where {T}
    # compute the magic square of odd orders
    A = zeros(T, n, n)
    i = 1
    j = div(n + 1, 2)
    for k = 1:n^2
        is = i
        js = j
        A[i, j] = k
        i = n - rem(n + 1 - i, n)
        j = rem(j, n) + 1
        if A[i, j] != 0
            i = rem(is, n) + 1
            j = js
        end
    end
    return A
end

"""
Magic Square Matrix
===================
The magic matrix is a matrix with integer entries such that
    the row elements, column elements, diagonal elements and
    anti-diagonal elements all add up to the same number.

*Input options:*

+ [type,] dim: the dimension of the matrix.

*Groups:* ["inverse"]
"""
function magic(::Type{T}, n::Integer) where {T}
    # Compute a magic square of order n
    # Learnt from Cleve Moler, Experiments with MATLAB, 2011
    if mod(n, 2) == 1
        # n is odd
        M = oddmagic(T, n)
    elseif mod(n, 4) == 0
        # n is doubly even
        a = floor.(Integer, mod.([1:n;], 4) / 2)
        B = broadcast(==, a', a)
        M = broadcast(+, T[1:n:n^2;]', T[0:n-1;])
        for i = 1:n, j = 1:n
            B[i, j] == 1 ? M[i, j] = n^2 + one(T) - M[i, j] : M[i, j]
        end
    else
        # n is singly even
        p = div(n, 2)
        M = oddmagic(T, p)
        M = [M M.+2*p^2; M.+3*p^2 M.+p^2]
        if n == 2
            return M
        end
        i = [1:p;]
        k = div(n - 2, 4)
        j = [[1:k;]; [(n-k+2):n;]]
        M[[i; i .+ p], j] = M[[i .+ p; i], j]
        i = k + 1
        j = [1, i]
        M[[i; i + p], j] = M[[i + p; i], j]
    end
    return M
end
magic(n::Integer) = magic(Float64, n)

"""
Grcar Matrix
============
The Grcar matrix is a Toeplitz matrix with sensitive
eigenvalues.

*Input options:*

+ [type,] dim, k: `dim` is the dimension of the matrix and
    `k` is the number of superdiagonals.

+ [type,] dim: the dimension of the matrix.

*Groups:* ["eigen"]

*References:*

**J. F. Grcar**, Operator coefficient methods
    for linear equations, Report SAND89-8691, Sandia National
    Laboratories, Albuquerque, New Mexico, 1989 (Appendix 2).
"""
function grcar(::Type{T}, n::Integer, k::Integer=3) where {T}
    # Compute grcar matrix
    G = tril(triu(ones(T, n, n)), min(k, n - 1)) - diagm(-1 => ones(T, n - 1))
    return G
end
grcar(args...) = grcar(Float64, args...)
grcar(::Type, args...) = throw(MethodError(grcar, Tuple(args)))

"""
Triw Matrix
===========
Upper triangular matrices discussed by Wilkinson and others.

*Input options:*

+ [type,] row_dim, col_dim, α, k: `row_dim` and `col_dim`
        are row and column dimension of the matrix. `α` is a
        scalar representing the entries on the superdiagonals.
        `k` is the number of superdiagonals.

+ [type,] dim: the dimension of the matrix.

*Groups:* ["inverse", "illcond"]

*References:*

**G. H. Golub and J. H. Wilkinson**, Ill-conditioned
eigensystems and the computation of the Jordan canonical form,
SIAM Review, 18(4), 1976, pp. 578-6
"""
function triw(::Type{T}, m::Integer, n::Integer, alpha, k::Integer) where {T}
    alpha = convert(T, alpha)
    A = tril(Matrix{T}(I, m, n) + alpha * triu(ones(T, m, n), 1), min(k, n - 1))
    return A
end
triw(::Type{T}, n::Integer) where {T} = triw(T, n, n, -1, n - 1)
triw(args...) = triw(Float64, args...)
triw(::Type, args...) = throw(MethodError(triw, Tuple(args)))

"""
Moler Matrix
============
The Moler matrix is a symmetric positive definite matrix.
It has one small eigenvalue.

*Input options:*

+ [type,] dim, alpha: `dim` is the dimension of the matrix,
        `alpha` is a scalar.

+ [type,] dim: `alpha = -1`.

*Groups:* ["inverse", "illcond", "symmetric", "posdef"]

*References:*

**J.C. Nash**, Compact Numerical Methods for Computers:
    Linear Algebra and Function Minimisation, second edition,
    Adam Hilger, Bristol, 1990 (Appendix 1).
"""
function moler(::Type{T}, n::Integer, alpha=-1.0) where {T}
    alpha = convert(T, alpha)
    M = triw(T, n, n, alpha, n - 1)' * triw(T, n, n, alpha, n - 1)
    return M
end
moler(args...) = moler(Float64, args...)
moler(::Type, args...) = throw(MethodError(moler, Tuple(args)))

"""
Pascal Matrix
=============
The Pascal matrix’s anti-diagonals form the Pascal’s triangle.

*Input options:*

+ [type,] dim: the dimension of the matrix.

*Groups:* ["inverse", "illcond", "symmetric", "posdef", "eigen"]

*References:*

**R. Brawer and M. Pirovino**, The linear algebra of
    the Pascal matrix, Linear Algebra and Appl., 174 (1992),
    pp. 13-23 (this paper gives a factorization of L = PASCAL(N,1)
               and a formula for the elements of L^k).

**N. J. Higham**, Accuracy and Stability of Numerical Algorithms,
second edition, Society for Industrial and Applied Mathematics, Philadelphia, PA,
USA, 2002; sec. 28.4.
"""
function pascal(::Type{T}, n::Integer) where {T}
    P = zeros(T, n, n)
    for j = 1:n
        for i = 1:n
            try
                P[i, j] = binomial(i + j - 2, j - 1)
            catch
                P[i, j] = binomial(BigInt(i + j - 2), j - 1)
            end
        end
    end
    return P
end
pascal(n::Integer) = pascal(Float64, n)

"""
Kahan Matrix
============
The Kahan matrix is an upper trapezoidal matrix, i.e., the
`(i,j)` element is equal to `0` if `i > j`. The useful range of
    `θ` is `0 < θ < π`. The diagonal is perturbed by
    `pert*eps()*diagm([n:-1:1;])`.

*Input options:*

+ [type,] rowdim, coldim, θ, pert: `rowdim` and `coldim` are the row and column
    dimensions of the matrix. `θ` and `pert` are scalars.

+ [type,] dim, θ, pert: `dim` is the dimension of the matrix.

+ [type,] dim: `θ = 1.2`, `pert = 25`.

*Groups:* ["inverse", "illcond"]

*References:*

**W. Kahan**, Numerical linear algebra, Canadian Math.
    Bulletin, 9 (1966), pp. 757-801.
"""
function kahan(::Type{T}, m::Integer, n::Integer, theta, pert) where {T}
    theta = convert(T, theta)
    pert = convert(T, pert)
    s = sin(theta)
    c = cos(theta)
    dim = min(m, n)
    A = zeros(T, m, n)
    for i = 1:m, j = 1:n
        i > dim ? A[i, j] = zero(T) :
        i > j ? A[i, j] = zero(T) :
        i == j ? A[i, j] = s^(i - 1) + pert * eps(T) * (m - i + 1) : A[i, j] = -c * s^(i - 1)
    end
    return A
end
kahan(::Type{T}, n::Integer, theta, pert) where {T} = kahan(T, n, n, theta, pert)
kahan(::Type{T}, n::Integer) where {T} = kahan(T, n, n, 1.2, 25.0)
kahan(args...) = kahan(Float64, args...)
kahan(::Type, args...) = throw(MethodError(kahan, Tuple(args)))

"""
Pei Matrix
==========
The Pei matrix is a symmetric matrix with known inversion.

*Input options:*

+ [type,] dim, α: `dim` is the dimension of the matrix.
    `α` is a scalar.

+ [type,] dim: the dimension of the matrix.

*Groups:* ["inverse", "illcond", "symmetric", "posdef"]

*References:*

**M. L. Pei**, A test matrix for inversion procedures,
    Comm. ACM, 5 (1962), p. 508.
"""
function pei(::Type{T}, n::Integer, alpha=1) where {T}
    alpha = convert(T, alpha)
    return alpha * Matrix{T}(I, n, n) + ones(T, n, n)
end
pei(args...) = pei(Float64, args...)
pei(::Type, args...) = throw(MethodError(pei, Tuple(args)))

"""
Vandermonde Matrix
==================
The inverse and determinat are known explicity.

*Input options:*

+ [type,] vec, col_dim: `vec` is a vector, `col_dim` is the number of columns.

+ [type,] vec: `col_dim = length(vec)`

+ [type,] dim: `vec = [1:dim;]`

*Groups:* ["inverse", "illcond"]

*References:*

**N. J. Higham**, Stability analysis of algorithms
    for solving confluent Vandermonde-like systems, SIAM J.
        Matrix Anal. Appl., 11 (1990), pp. 23-41.
"""
function vand(::Type{T}, p::Vector, n::Integer) where {T}
    # n: number of rows
    # p: a vector
    m = length(p)
    V = Array{T,2}(undef, m, n)
    for j = 1:m
        @inbounds V[j, 1] = 1
    end
    for i = 2:n
        for j = 1:m
            @inbounds V[j, i] = p[j] * V[j, i-1]
        end
    end
    return V
end
vand(::Type{T}, n::Integer) where {T} = vand(T, [1:n;], n)
vand(::Type{T}, p::Vector) where {T} = vand(T, p, length(p))
vand(args...) = vand(Float64, args...)
vand(::Type, args...) = throw(MethodError(vand, Tuple(args)))

"""
Involutory Matrix
=================
An involutory matrix is a matrix that is its own inverse.

*Input options:*

+ [type,] dim: `dim` is the dimension of the matrix.

*Groups:* ["inverse", "illcond", "eigen"]

*References:*

**A. S. Householder and J. A. Carpenter**, The
        singular values of involutory and of idempotent matrices,
        Numer. Math. 5 (1963), pp. 234-237.
"""
function invol(::Type{T}, n::Integer) where {T}
    A = hilb(T, n)
    d = -n
    A[:, 1] = d * A[:, 1]
    for i = 1:n-1
        d = -(n + i) * (n - i) * d / (i * i)
        A[i+1, :] = d * A[i+1, :]
    end
    return A
end
invol(n::Integer) = invol(Float64, n)

"""
Chebyshev Spectral Differentiation Matrix
=========================================
If `k = 0`,the generated matrix is nilpotent and a vector with
        all one entries is a null vector. If `k = 1`, the generated
        matrix is nonsingular and well-conditioned. Its eigenvalues
        have negative real parts.

*Input options:*

+ [type,] dim, k: `dim` is the dimension of the matrix and
        `k = 0 or 1`.

+ [type,] dim: `k=0`.

*Groups:* ["eigen"]

*References:*

**L. N. Trefethen and M. R. Trummer**, An instability
        phenomenon in spectral methods, SIAM J. Numer. Anal., 24 (1987), pp. 1008-1023.
"""
function chebspec(::Type{T}, n::Integer, k::Integer=0) where {T}
    # k = 0 or 1
    if k == 1
        n = n + 1
    end
    c = ones(T, n)
    c[1] = 2
    c[2:n-1] .= 1
    c[n] = 2
    x = ones(T, n)

    A = zeros(T, n, n)
    # Compute the chebyshev points
    for i = 1:n
        x[i] = cos(pi * (i - 1) / (n - 1))
    end

    for j = 1:n, i = 1:n
        i != j ? A[i, j] = (-1)^(i + j) * c[i] / (c[j] * (x[i] - x[j])) :
        i == 1 ? A[i, i] = (2 * (n - 1)^2 + 1) / 6 :
        i == n ? A[i, i] = -(2 * (n - 1)^2 + 1) / 6 :
        A[i, i] = -0.5 * x[i] / (1 - x[i]^2)
    end
    if k == 1
        A = A[2:n, 2:n]
    end
    return A
end
chebspec(args...) = chebspec(Float64, args...)
chebspec(::Type, args...) = throw(MethodError(chebspec, Tuple(args)))

"""
Lotkin Matrix
=============
The Lotkin matrix is the Hilbert matrix with its first row
        altered to all ones. It is unsymmetric, illcond and
        has many negative eigenvalues of small magnitude.

*Input options:*

+ [type,] dim: `dim` is the dimension of the matrix.

*Groups:* ["inverse", "illcond", "eigen"]

*References:*

**M. Lotkin**, A set of test matrices, MTAC, 9 (1955), pp. 153-161.
"""
function lotkin(::Type{T}, n::Integer) where {T}
    A = hilb(T, n)
    A[1, :] = ones(T, n)'
    return A
end
lotkin(n::Integer) = lotkin(Float64, n)

"""
Clement Matrix
==============
The Clement matrix is a tridiagonal matrix with zero
        diagonal entries. If k = 1, the matrix is symmetric.

*Input options:*

+ [type,] dim, k: `dim` is the dimension of the matrix.
        If `k = 0`, the matrix is of type `Tridiagonal`.
        If `k = 1`, the matrix is of type `SymTridiagonal`.

+ [type,] dim: `k = 0`.

*Groups:* ["inverse", "symmetric", "eigen"]

*References:*

**P. A. Clement**, A class of triple-diagonal
        matrices for test purposes, SIAM Review, 1 (1959), pp. 50-52.
"""
function clement(::Type{T}, n::Integer, k::Integer=0) where {T}
    # construct Tridiagonal matrix
    # n is the dimension of the matrix
    # k = 0 or 1
    if n == 1 # handle the 1-d case.
        return zeros(T, 1, 1)
    end
    n = n - 1
    x = T[n:-1:1;]
    z = T[1:n;]
    if k == 0
        A = Tridiagonal(x, zeros(T, n + 1), z)
    else
        y = sqrt.(x .* z)
        A = SymTridiagonal(zeros(T, n + 1), y)
    end
    return A
end
clement(args...) = clement(Float64, args...)
clement(::Type, args...) = throw(MethodError(clement, Tuple(args)))

"""
Fiedler Matrix
==============
The Fiedler matrix is symmetric matrix with a dominant
      positive eigenvalue and all the other eigenvalues are negative.

*Input options:*

+ [type,] vec: a vector.

+ [type,] dim: `dim` is the dimension of the matrix. `vec=[1:dim;]`.

*Groups: *["inverse", "symmetric", "eigen"]

*References:*

**G. Szego**, Solution to problem 3705, Amer. Math.
            Monthly, 43 (1936), pp. 246-259.

**J. Todd**, Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
            Birkhauser, Basel, and Academic Press, New York, 1977, p. 159.
"""
function fiedler(::Type{T}, v::Vector) where {T}
    n = length(v)
    v = transpose(v[:])
    A = ones(T, n) * v
    A = abs.(A - transpose(A)) # nonconjugate transpose
end
fiedler(::Type{T}, n::Integer) where {T} = fiedler(T, [1:n;])
fiedler(args...) = fiedler(Float64, args...)
fiedler(::Type, args...) = throw(MethodError(fiedler, Tuple(args)))

"""
MIN[I,J] Matrix
===============
A matrix with `(i,j)` entry `min(i,j)`. It is a symmetric positive
     definite matrix. The eigenvalues and eigenvectors are known
     explicitly. Its inverse is tridiagonal.

*Input options:*

+ [type,] dim: the dimension of the matrix.

*Groups:* ["inverse", "symmetric", "posdef", "eigen"]

*References:*

**J. Fortiana and C. M. Cuadras**, A family of matrices,
            the discretized Brownian bridge, and distance-based regression,
            Linear Algebra Appl., 264 (1997), 173-188.  (For the eigensystem of A.)
"""
function minij(::Type{T}, n::Integer) where {T}
    A = zeros(T, n, n)
    [A[i, j] = min(i, j) for i = 1:n, j = 1:n]
    return A
end
minij(n::Integer) = minij(Float64, n)

"""
Binomial Matrix
===============
The matrix is a multiple of an involutory matrix.

*Input options:*

+ [type,] dim: the dimension of the matrix.
"""
function binomialm(::Type{T}, n::Integer) where {T}
    # Mulitiple of involutory matrix
    L = Array{T,2}(undef, n, n)
    D = Diagonal((-2) .^ [0:n-1;])
    [L[i, j] = binomial(i - 1, j - 1) for i = 1:n, j = 1:n]
    U = L[n:-1:1, n:-1:1]
    return L * D * U
end
binomialm(n::Integer) = binomialm(Float64, n)

"""
Tridiagonal Matrix
==================
Construct a tridigonal matrix of type `Tridiagonal`.

*Input options:*

+ [type,] v1, v2, v3: `v1` and `v3` are vectors of subdiagonal
            and superdiagonal elements, respectively, and `v2` is a vector
            of diagonal elements.

+ [type,] dim, x, y, z: `dim` is the dimension of the matrix,
            `x`, `y`, `z` are scalars. `x` and `z` are the subdiagonal and
            superdiagonal elements, respectively, and `y` is the diagonal
            elements.

+ [type,] dim: `x = -1, y = 2, z = -1`. This matrix is also
            known as the second difference matrix.

*Groups:* ["inverse", "illcond", "posdef", "eigen"]

*References:*

**J. Todd**, Basic Numerical Mathematics, Vol. 2:
            Numerical Algebra, Birkhauser, Basel, and Academic Press,
            New York, 1977, p. 155.
"""
function tridiag(::Type{T}, x::AbstractVector, y::AbstractVector,
    z::AbstractVector) where {T}
    x = map((i) -> convert(T, i), x)
    y = map((i) -> convert(T, i), y)
    z = map((i) -> convert(T, i), z)
    return Tridiagonal(x, y, z)
end
# Toeplitz tridiagonal matrix
tridiag(::Type{T}, n::Integer, x::Integer, y::Integer, z::Integer) where {T} =
    n == 1 ? y * ones(T, 1, 1) :
    tridiag(T, x * ones(T, n - 1), y * ones(T, n), z * ones(T, n - 1))
tridiag(::Type{T}, n::Integer) where {T} = tridiag(T, n, -1, 2, -1)
tridiag(args...) = tridiag(Float64, args...)
tridiag(::Type, args...) = throw(MethodError(tridiag, Tuple(args)))

"""
Lehmer Matrix
=============
The Lehmer matrix is a symmetric positive definite matrix.
            It is totally nonnegative. The inverse is tridiagonal and
            explicitly known

*Input options:*

+ [type,] dim: the dimension of the matrix.

*Groups:* ["inverse", "symmetric", "posdef"]

*References:*

**M. Newman and J. Todd**, The evaluation of
            matrix inversion programs, J. Soc. Indust. Appl. Math.,
            6 (1958), pp. 466-476.
            Solutions to problem E710 (proposed by D.H. Lehmer): The inverse
            of a matrix, Amer. Math. Monthly, 53 (1946), pp. 534-535.
"""
function lehmer(::Type{T}, n::Integer) where {T}
    A = Array{T,2}(undef, n, n)
    [A[i, j] = min(i, j) / max(i, j) for i = 1:n, j = 1:n]
    return A
end
lehmer(n::Integer) = lehmer(Float64, n)
#=
function lehmer(::Type{T}, n::Integer) where T
    A = Array{T,2}(n, n)
    [A[i,j] = min(i,j) / max(i,j) for i = 1:n, j = 1:n]
    return A
end
lehmer(n::Integer) = lehmer(Float64, n)
=#
"""
Parter Matrix
=============
The Parter matrix is a Toeplitz and Cauchy matrix
            with singular values near `π`.

*Input options:*

+ [type,] dim: the dimension of the matrix.

*Groups:* ["eigen"]

*References:*

The MathWorks Newsletter, Volume 1, Issue 1,
            March 1986, page 2. S. V. Parter, On the distribution of the
            singular values of Toeplitz matrices, Linear Algebra and
            Appl., 80 (1986), pp. 115-130.
"""
function parter(::Type{T}, n::Integer) where {T}
    A = Array{T,2}(undef, n, n)
    [A[i, j] = one(T) / (i - j + 0.5) for i = 1:n, j = 1:n]
    return A
end
parter(n::Integer) = parter(Float64, n)

"""
Chow Matrix
===========
The Chow matrix is a singular Toeplitz lower Hessenberg matrix.

*Input options:*

+ [type,] dim, alpha, delta: `dim` is dimension of the matrix.
            `alpha`, `delta` are scalars such that `A[i,i] = alpha + delta` and
            `A[i,j] = alpha^(i + 1 -j)` for `j + 1 <= i`.

+ [type,] dim: `alpha = 1`, `delta = 0`.

*Groups:* ["eigen"]

*References:*

**T. S. Chow**, A class of Hessenberg matrices with known
                eigenvalues and inverses, SIAM Review, 11 (1969), pp. 391-395.
"""
function chow(::Type{T}, n::Integer, alpha, delta) where {T}
    A = zeros(T, n, n)
    alpha = convert(T, alpha)
    delta = convert(T, delta)
    for i = 1:n, j = 1:n
        if i == j - 1
            A[i, j] = one(T)
        elseif i == j
            A[i, j] = alpha + delta
        elseif j + 1 <= i
            A[i, j] = alpha^(i + 1 - j)
        end
    end
    return A
end
chow(::Type{T}, n::Integer) where {T} = chow(T, n, 1, 0)
chow(args...) = chow(Float64, args...)
chow(::Type, args...) = throw(MethodError(chow, Tuple(args)))

#
# newsign: newsign(0) = 1
#
function newsign(x)
    x == 0 ? y = 1 : y = sign(x)
    return y
end

"""
Random Correlation Matrix
=========================
A random correlation matrix is a symmetric positive
     semidefinite matrix with 1s on the diagonal.

*Input options:*

+ [type,] dim [; rng]: the dimension of the matrix. `rng` is an `AbstractRNG`

*Groups:* ["symmetric", 'pos-semidef', "random"]
"""
function randcorr(::Type{T}, n::Integer; rng=Random.default_rng()) where {T}
    rng::AbstractRNG
    x = rand(rng, T, n) # x is the vector of random eigenvalues from a uniform distribution.
    x = n * x / sum(x) # x has nonnegtive elements.
    A = diagm(0 => x)
    F = qr(randn(rng, n, n))
    Q = F.Q * diagm(0 => sign.(diag(F.R))) # form a random orthogonal matrix.
    copyto!(A, Q * A * Q')

    a = diag(A)
    l = findall(a .< 1)
    g = findall(a .> 1)

    # Apply Given rotation to set A[i,i] = 1
    while length(l) > 0 && length(g) > 0
        k = ceil(Integer, rand(rng) * length(l))
        h = ceil(Integer, rand(rng) * length(g))
        i = l[k]
        j = g[h]
        if i > j
            i, j = j, i
        end
        alpha = sqrt(A[i, j]^2 - (a[i] - 1) * (a[j] - 1))
        # take sign to avoid cancellation.
        t = (A[i, j] + newsign(A[i, j]) * alpha) / (a[j] - 1)
        c = 1 / sqrt(1 + t^2)
        s = t * c

        A[:, [i, j]] = A[:, [i, j]] * [c s; -s c]
        A[[i, j], :] = [c -s; s c] * A[[i, j], :]

        A[i, i] = 1
        a = diag(A)
        l = findall(a .< 1)
        g = findall(a .> 1)
    end
    [A[i, i] = 1 for i = 1:n]
    return (A + A') / 2

end
randcorr(args...; kw...) = randcorr(Float64, args...; kw...)
randcorr(::Type, args...; kw...) = throw(MethodError(randcorr, Tuple(args)))

"""
Poisson Matrix
==============
A block tridiagonal matrix from Poisson’s equation.
     This matrix is sparse, symmetric positive definite and
     has known eigenvalues. The result is of type `SparseMatrixCSC`.

*Input options:*

+ [type,] dim: the dimension of the matirx is `dim^2`.

*Groups:* ["inverse", "symmetric", "posdef", "eigen", "sparse"]

*References:*

**G. H. Golub and C. F. Van Loan**, Matrix Computations,
          second edition, Johns Hopkins University Press, Baltimore,
          Maryland, 1989 (Section 4.5.4).
"""
function poisson(::Type{T}, n::Integer) where {T}
    S = Array(tridiag(T, n))
    A = sparse(T(1)I, n, n)
    return kron(A, S) + kron(S, A)
end
poisson(n::Integer) = poisson(Float64, n)

"""
Toeplitz Matrix
===============
A Toeplitz matrix is a matrix in which each descending
       diagonal from left to right is constant.

*Input options:*

+ [type,] vc, vr: `vc` and `vr` are the first column and row of the matrix.

+ [type,] v: symmatric case, i.e., `vc = vr = v`.

+ [type,] dim: `dim` is the dimension of the matrix. `v = [1:dim;]` is the first
                row and column vector.
"""
function toeplitz(::Type{T}, vc::Vector, vr::Vector) where {T}
    n = length(vc)
    length(vr) == n || throw(DimensionMismatch(""))
    vc[1] == vr[1] || error("The first element of the vectors must be the same.")
    A = Array{T,2}(undef, n, n)
    [i >= j ? A[i, j] = vc[i-j+1] : A[i, j] = vr[j-i+1] for i = 1:n, j = 1:n]
    A
end
toeplitz(::Type{T}, v::Vector) where {T} = toeplitz(T, v, v)
toeplitz(::Type{T}, n::Integer) where {T} = toeplitz(T, [1:n;])
toeplitz(args...) = toeplitz(Float64, args...)
toeplitz(::Type, args...) = throw(MethodError(toeplitz, Tuple(args)))

"""
Hankel Matrix
=============
A Hankel matrix is a matrix that is symmetric and constant
                across the anti-diagonals.

*Input options:*

+ [type,] vc, vr: `vc` and `vc` are the first column and last row of the
       matrix. If the last element of `vc` differs from the first element
                of `vr`, the last element of `rc` prevails.

+ [type,] v: `vc = vr = v`.

+ [type,] dim: `dim` is the dimension of the matrix. `v = [1:dim;]`.
"""
function hankel(::Type{T}, vc::Vector, vr::Vector) where {T}
    p = [vc; vr[2:end]]
    m = length(vc)
    n = length(vr)
    H = Array{T,2}(undef, m, n)
    [H[i, j] = p[i+j-1] for i = 1:m, j = 1:n]
    H
end
hankel(::Type{T}, v::Vector) where {T} = hankel(T, v, v)
hankel(::Type{T}, n::Integer) where {T} = hankel(T, [1:n;])
hankel(args...) = hankel(Float64, args...)
hankel(::Type, args...) = throw(MethodError(hankel, Tuple(args)))

"""
Prolate Matrix
==============
A prolate matrix is a symmetirc, ill-conditioned Toeplitz matrix.

*Input options:*

+ [type,] dim, w: `dim` is the dimension of the matrix. `w` is a real scalar.

+ [type,] dim: the case when `w = 0.25`.

*References:*

**J. M. Varah**. The Prolate Matrix. Linear Algebra and Appl.
             187:267--278, 1993.
"""
function prolate(::Type{T}, n::Integer, w::Real) where {T}
    v = Array{T,1}(undef, n)
    v[1] = 2 * w
    [v[i] = sin(2 * pi * w * i) / pi * i for i = 2:n]
    return toeplitz(T, v)
end
prolate(::Type{T}, n::Integer) where {T} = prolate(T, n, 0.25)
prolate(args...) = prolate(Float64, args...)
prolate(::Type, args...) = throw(MethodError(prolate, Tuple(args)))

"""
Neumann Matrix
==============
A singular matrix from the discrete Neumann problem.
       The matrix is sparse and the null space is formed by a vector of ones

*Input options:*

+ [type,] dim: the dimension of the matrix is `dim^2`.

*Groups:* ["eigen", "sparse"]

*References:*

**R. J. Plemmons**, Regular splittings and the
          discrete Neumann problem, Numer. Math., 25 (1976), pp. 153-161.
"""
function neumann(::Type{T}, n::Integer) where {T}
    if n == 1
        return 4 * ones(T, 1, 1) #handle 1-d case.
    end
    S = Matrix(tridiag(T, n))
    S[1, 2] = -2
    S[n, n-1] = -2
    A = sparse(T(1)I, n, n)
    return kron(S, A) + kron(A, S)
end
neumann(n::Integer) = neumann(Float64, n)

#
# Sylvester's orthogonal matrix
# See Rosser matrix References 2.
#
# for a = d = 2, b = c = 1, P_block' * P_block = 10 * Identity
#
P_block(::Type{T}, a, b, c, d) where {T} =
    reshape(T[a, b, c, d, b, -a, -d, c, c, d, -a, -b, d, -c, b, -a], 4, 4)

"""
Rosser Matrix
=============
The Rosser matrix’s eigenvalues are very close together
        so it is a challenging matrix for many eigenvalue algorithms.

*Input options:*

+ [type,] dim, a, b: `dim` is the dimension of the matrix.
            `dim` must be a power of 2.
            `a` and `b` are scalars. For `dim = 8, a = 2` and `b = 1`, the generated
            matrix is the test matrix used by Rosser.

+ [type,] dim: `a = rand(1:5), b = rand(1:5)`.

*Groups:* ["eigen", "illcond", "random"]

*References:*

**J. B. Rosser, C. Lanczos, M. R. Hestenes, W. Karush**,
            Separation of close eigenvalues of a real symmetric matrix,
            Journal of Research of the National Bureau of Standards, v(47)
            (1951)
"""
function rosser(::Type{T}, n::Integer, a, b) where {T}
    if n < 1
        lgn = 0
    else
        lgn = round(Integer, log2(n))
    end
    2^lgn != n && throw(ArgumentError("n must be positive integer and a power of 2."))

    if n == 1 # handle 1-d case
        return 611 * ones(T, 1, 1)
    end

    if n == 2
        #eigenvalues are 500, 510
        B = T[101 1; 1 101]
        P = T[2 1; 1 -2]
        A = P' * B * P
    elseif n == 4
        # eigenvalues are 0.1, 1019.9, 1020, 1020 for a = 2 and b = 1
        B = zeros(T, n, n)
        B[1, 1], B[1, 4], B[4, 1], B[4, 4] = 101, 1, 1, 101
        B[2, 2], B[2, 3], B[3, 2], B[3, 3] = 1, 10, 10, 101
        P = P_block(T, a, b, b, a)
        A = P' * B * P
    elseif n == 8
        # eigenvalues are 1020, 1020, 1000, 1000, 0.098, 0, -1020
        B = zeros(T, n, n)
        B[1, 1], B[6, 1], B[2, 2], B[8, 2] = 102, 1, 101, 1
        B[3, 3], B[7, 3] = 98, 14
        B[4, 4], B[5, 4], B[4, 5], B[5, 5] = 1, 10, 10, 101
        B[1, 6], B[6, 6], B[3, 7], B[7, 7], B[2, 8], B[8, 8] = 1, -102, 14, 2, 1, 101
        P = [P_block(T, a, b, b, a)' zeros(T, 4, 4); zeros(T, 4, 4) P_block(T, b, -b, -a, a)]
        A = P' * B * P
    else
        lgn = lgn - 2
        halfn = round(Integer, n / 2)
        # using Sylvester's method
        P = P_block(T, a, b, b, a)
        m = 4
        for i in 1:lgn
            P = [P zeros(T, m, m); zeros(T, m, m) P]
            m = m * 2
        end
        # mix 4 2-by-2 matrices (with close eigenvalues) into a large nxn matrix.
        B_list = T[102, 1, 1, -102, 101, 1, 1, 101, 1, 10, 10, 101, 98, 14, 14, 2]
        B = zeros(T, n, n)
        j, k = 1, 5
        for i in 1:(halfn+1)
            indexend = halfn - 1 + i
            list_start = k
            list_end = k + 3

            if list_start > 16 || list_end > 16
                k = 1
                list_start = 1
                list_end = 4
            end
            B[j, j], B[j, indexend], B[indexend, j], B[indexend, indexend] = B_list[list_start:list_end]
            j = j + 1
            k = k + 4
        end
        A = P' * B * P
    end
    symmetrize!(A)
    return A
end
function rosser(::Type{T}, n::Integer; rng=Random.default_rng()) where {T}
    rosser(T, n, rand(rng, 1:5), rand(rng, 1:5))
end
rosser(args...; kw...) = rosser(Float64, args...; kw...)
rosser(::Type, args...; kw...) = throw(MethodError(rosser, Tuple(args)))

"""
    symmetrize!(A::Matrix)

    Force `A` to be symmetric real/hermitian
"""
function symmetrize!(A::AbstractMatrix)
    m, n = size(A)
    for j = 1:n
        for i = 1:j-1
            aij, aji = A[i, j], conj(A[j, i])
            if aij != aji
                aij = aij + (aji - aij) / 2
                A[i, j] = aij
                A[j, i] = conj(aij)
            end
        end
        A[j,j] = real(A[j,j])
    end
    A
end

"""
Matrix with Application in Sampling Theory
==========================================
A nonsymmetric matrix with eigenvalues 0, 1, 2, ... n-1.

*Input options:*

+ [type,] vec: `vec` is a vector with no repeated elements.

+ [type,] dim: `dim` is the dimension of the matrix.
            `vec = [1:dim;]/dim`.

*Groups:* ["eigen"]

*References:*

**L. Bondesson and I. Traat**, A nonsymmetric matrix
            with integer eigenvalues, linear and multilinear algebra, 55(3)
            (2007), pp. 239-247
"""
function sampling(::Type{T}, x::Vector) where {T}
    n = length(x)
    A = zeros(T, n, n)
    for j = 1:n, i = 1:n
        if i != j
            A[i, j] = x[i] / (x[i] - x[j])
        end
    end
    d = sum(A, dims=2)
    A = A + diagm(0 => d[:])
    return A
end
#
# special probability case
# see:
#   L. Bondesson and I. Traat, A Nonsymmetric Matrix with Integer
#   Eigenvalues, Linear and Multilinear Algebra, 55(3)(2007), pp. 239-247.
#
function sampling(::Type{T}, n::Integer) where {T}
    p = T[1:n;] / n
    return sampling(T, p)
end
sampling(args...) = sampling(Float64, args...)
sampling(::Type, args...) = throw(MethodError(sampling, Tuple(args)))

"""
Wilkinson Matrix
================
The Wilkinson matrix is a symmetric tridiagonal matrix with pairs
of nearly equal eigenvalues. The most frequently used case
is `matrixdepot("wilkinson", 21)`. The result is of type `Tridiagonal`.

*Input options:*

+ [type,] dim: the dimension of the matrix.

*Groups:* ["symmetric", "eigen"]

*References:*

**J. H. Wilkinson**, Error analysis of direct methods
of matrix inversion, J. Assoc. Comput. Mach., 8 (1961),  pp. 281-330.
"""
function wilkinson(::Type{T}, n::Integer) where {T}
    if n == 1 # handle 1-d case
        return ones(T, 1, 1)
    end
    m = (n - 1) / 2
    A = Tridiagonal(ones(T, n - 1), abs.(T[-m:m;]), ones(T, n - 1))
    return A
end
wilkinson(n::Integer) = wilkinson(Float64, n)

using Random

"""
Random Matrix with Element -1, 0, 1
===================================

*Input options:*

+ [type,] row_dim, col_dim, k[; rng]: `row_dim` and `col_dim` are row and column dimensions,
   `rng` is an `AbstractRNG` with default
   `k = 1`: entries are 0 or 1.
   `k = 2`: entries are -1 or 1.
   `k = 3`: entries are -1, 0 or 1.

+ [type,] dim, k: `row_dim = col_dim = dim`.

+ [type,] dim: `k = 1`.

*Groups:* ["random"]
"""
function rando(::Type{T}, m::Integer, n::Integer, k::Integer;
    rng=Random.default_rng()) where {T}

    rng::AbstractRNG
    A = Array{T,2}(undef, m, n)
    if k == 1
        copyto!(A, floor.(rand(rng, m, n) .+ 0.5))
    elseif k == 2
        copyto!(A, 2 * floor.(rand(rng, m, n) .+ 0.5) .- one(T))
    elseif k == 3
        copyto!(A, round.(3 * rand(rng, m, n) .- 1.5))
    else
        throw(ArgumentError("invalid k value."))
    end
    return A
end
rando(::Type{T}, n::Integer, k::Integer; kw...) where {T} = rando(T, n, n, k; kw...)
rando(::Type{T}, n::Integer; kw...) where {T} = rando(T, n, n, 1; kw...)
rando(args...; kw...) = rando(Float64, args...; kw...)
rando(::Type, args...; kw...) = throw(MethodError(rando, Tuple(args)))

#
# Pre-multiply by random orthogonal matrix
#
function qmult!(rng::AbstractRNG, A::Matrix{T}) where {T}
    n, m = size(A)

    d = zeros(T, n)
    for k = n-1:-1:1

        # generate random Householder transformation
        x = randn(rng, n - k + 1)
        s = norm(x)
        sgn = sign(x[1]) + (x[1] == 0)
        s = sgn * s
        d[k] = -sgn
        x[1] = x[1] + s
        beta = s * x[1]

        # apply the transformation to A
        y = x' * A[k:n, :]
        A[k:n, :] = A[k:n, :] - x * (y / beta)
    end

    # tidy up signs
    for i = 1:n-1
        A[i, :] = d[i] * A[i, :]
    end
    A[n, :] = A[n, :] * sign(randn(rng))
    return A
end

"""
Random Matrix with Pre-assigned Singular Values
===============================================
*Input options:*

+ [type,] row_dim, col_dim, sigma::Vector [; rng]
  `sigma` is the vector of singular values
  `rng`   is an `AbstractRNG` with default

+ [type,] row_dim, col_dim, κ, mode: `row_dim` and `col_dim`
    are the row and column dimensions.
  `κ` is the condition number of the matrix.
  `mode = 1`: one large singular value.
  `mode = 2`: one small singular value.
  `mode = 3`: geometrically distributed singular values.
  `mode = 4`: arithmetrically distributed singular values.
  `mode = 5`: random singular values with  unif. dist. logarithm.

+ [type,] dim, κ, mode: `row_dim = col_dim = dim`.

+ [type,] dim, κ: `mode = 3`.

+ [type,] dim: `κ = sqrt(1/eps())`, `mode = 3`.

*Groups:* ["illcond", "random"]

*References:*

**N. J. Higham**, Accuracy and Stability of Numerical
Algorithms, second edition, Society for Industrial and Applied Mathematics,
Philadelphia, PA, USA, 2002; sec. 28.3.
"""
function randsvd(::Type{T}, m::Integer, n::Integer, κ::AbstractFloat, mode::Integer=3;
    rng=Random.default_rng()) where {T}

    sigma = make_sigma(rng, min(m, n), T(κ), mode)
    randsvd(T, m, n, sigma; rng)
end

function randsvd(::Type{T}, m::Integer, n::Integer, sigma::AbstractVector;
    rng=Random.default_rng()) where {T}

    rng::AbstractRNG
    A = zeros(T, m, n)
    p = length(sigma)
    for i in 1:min(length(sigma), m, n)
        A[i, i] = sigma[i]
    end
    A = qmult!(rng, copy(A'))
    A = qmult!(rng, copy(A'))

    return A
end
function randsvd(::Type{T}, n::Integer, κ::AbstractFloat, mode::Integer=3; kw...) where {T}
    randsvd(T, n, n, κ, mode; kw...)
end
randsvd(::Type{T}, n::Integer; kw...) where {T} = randsvd(T, n, sqrt(1 / eps(T)); kw...)
randsvd(args...; kw...) = randsvd(Float64, args...; kw...)
randsvd(::Type, args...; kw...) = throw(MethodError(randsvd, Tuple(args)))

function make_sigma(rng::AbstractRNG, p::Integer, κ::T, mode::Integer) where {T}

    abs(κ) >= 1 || throw(ArgumentError("Condition number must be at least 1."))

    if p == 1 # handle 1-d case
        return ones(T, p) * κ
    end

    sigma = if mode == 3
        factor = κ^(-1 / (p - 1))
        factor .^ [0:p-1;]
    elseif mode == 4
        ones(T, p) - T[0:p-1;] / (p - 1) * (1 - 1 / κ)
    elseif mode == 5
        exp.(-rand(rng, p) * log(κ))
    elseif mode == 2
        sigmax = ones(T, p)
        sigmax[p] = one(T) / κ
        sigmax
    elseif mode == 1
        sigmax = ones(p) ./ κ
        sigmax[1] = one(T)
        sigmax
    else
        throw(ArgumentError("invalid mode value."))
    end
    return sigma
end

"""
Random Orthogonal Upper Hessenberg Matrix
=========================================
The matrix is constructed via a product of Givens rotations.

*Input options:*

+ [type,] dim: the dimension of the matrix.

*Groups:* ["random"]

*References:*

**W. B. Gragg**, The QR algorithm for unitary
    Hessenberg matrices, J. Comp. Appl. Math., 16 (1986), pp. 1-8.
"""
function rohess(::Type{T}, n::Integer; rng=Random.default_rng()) where {T}
    rng::AbstractRNG
    x = rand(rng, n - 1) * 2 * pi
    H = Matrix{T}(I, n, n)
    H[n, n] = sign(randn(rng))
    for i = n:-1:2
        theta = x[i-1]
        c = convert(T, cos(theta))
        s = convert(T, sin(theta))
        H[[i - 1; i], :] = [c * H[i-1, :][:]' + s * H[i, :][:]'; -s * H[i-1, :][:]' + c * H[i, :][:]']
    end
    return H
end
rohess(n::Integer; kw...) = rohess(Float64, n, kw...)

"""
Kac-Murdock-Szego Toeplitz matrix
=================================

*Input options:*

+ [type,] dim, rho: `dim` is the dimension of the matrix, `rho` is a
    scalar such that `A[i,j] = rho^(abs(i-j))`.

+ [type,] dim: `rho = 0.5`.

*Groups:* ["inverse", "illcond", "symmetric", "posdef"]

*References:*

**W. F. Trench**, Numerical solution of the eigenvalue
    problem for Hermitian Toeplitz matrices, SIAM J. Matrix Analysis
    and Appl., 10 (1989), pp. 135-146 (and see the references therein).
"""
function kms(::Type{T}, n::Integer, rho::Number) where {T}
    A = typeof(rho) <: Complex ? Array{typeof(rho)}(undef, n, n) : Array{T,2}(undef, n, n)
    [A[i, j] = rho^(abs(i - j)) for i = 1:n, j = 1:n]
    if typeof(rho) <: Complex
        A = conj(tril(A, -1)) + triu(A)
    end
    return A
end
kms(::Type{T}, n::Integer) where {T} = kms(T, n, convert(T, 0.5))
kms(args...) = kms(Float64, args...)
kms(::Type, args...) = throw(MethodError(kms, Tuple(args)))

"""
Wathen Matrix
=============
Wathen Matrix is a sparse, symmetric positive, random matrix
arose from the finite element method. The generated matrix is
the consistent mass matrix for a regular nx-by-ny grid of
8-nodes.

*Input options:*

+ [type,] nx, ny: the dimension of the matrix is equal to
    `3 * nx * ny + 2 * nx * ny + 1`.

+ [type,] n: `nx = ny = n`.

*Groups:* ["symmetric", "posdef", "eigen", "random", "sparse"]

*References:*

**A. J. Wathen**, Realistic eigenvalue bounds for
    the Galerkin mass matrix, IMA J. Numer. Anal., 7 (1987),
    pp. 449-457.
"""
function wathen(::Type{T}, nx::Integer, ny::Integer; rng=Random.default_rng()) where {T}

    rng::AbstractRNG
    e1 = T[6 -6 2 -8; -6 32 -6 20; 2 -6 6 -6; -8 20 -6 32]
    e2 = T[3 -8 2 -6; -8 16 -8 20; 2 -8 3 -8; -6 20 -8 16]
    e3 = [e1 e2; e2' e1] / 45
    n = 3 * nx * ny + 2 * nx + 2 * ny + 1
    ntriplets = nx * ny * 64
    Irow = zeros(Int, ntriplets)
    Jrow = zeros(Int, ntriplets)
    Xrow = zeros(T, ntriplets)
    ntriplets = 0
    rho = 100 * rand(rng, nx, ny)
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

            em = convert(T, rho[i, j]) * e3

            for krow = 1:8
                for kcol = 1:8
                    ntriplets += 1
                    Irow[ntriplets] = node[krow]
                    Jrow[ntriplets] = node[kcol]
                    Xrow[ntriplets] = em[krow, kcol]
                end
            end

        end
    end
    return sparse(Irow, Jrow, Xrow, n, n)
end
wathen(::Type{T}, n::Integer; kw...) where {T} = wathen(T, n, n; kw...)
wathen(args...; kw...) = wathen(Float64, args...; kw...)
wathen(::Type, args...; kw...) = throw(MethodError(wathen, Tuple(args)))

"""
Golub Matrix
============
Golub matrix is the product of two random unit lower and upper
    triangular matrices respectively. LU factorization without pivoting
    fails to reveal that such matrices are badly conditioned.

*Input options:*

+ [type,] dim: the dimension of the matrix.

*References:*

**D. Viswanath and N. Trefethen**. Condition Numbers of
    Random Triangular Matrices, SIAM J. Matrix Anal. Appl. 19, 564-581,
    1998.
"""
function golub(::Type{T}, n::Integer; rng=Random.default_rng()) where {T}
    rng::AbstractRNG
    s = 10
    L = Array{T,2}(undef, n, n)
    U = Array{T,2}(undef, n, n)
    if T <: Integer
        [L[i, j] = round_matlab(T, s * randn(rng)) for j = 1:n, i = 1:n]
        [U[i, j] = round_matlab(T, s * randn(rng)) for j = 1:n, i = 1:n]
    else
        [L[i, j] = s * randn(rng) for j = 1:n, i = 1:n]
        [U[i, j] = s * randn(rng) for j = 1:n, i = 1:n]
    end
    L = tril(L, -1) + Matrix{T}(I, n, n)
    U = triu(U, 1) + Matrix{T}(I, n, n)
    return L * U
end
golub(n::Integer; kw...) = golub(Float64, n; kw...)

"""
Companion Matrix
================
The companion matrix to a monic polynomial
    `a(x) = a_0 + a_1x + ... + a_{n-1}x^{n-1} + x^n`
    is the n-by-n matrix with ones on the subdiagonal and
    the last column given by the coefficients of `a(x)`.

*Input options:*

+ [type,] vec: `vec` is a vector of coefficients.

+ [type,] dim: `vec = [1:dim;]`. `dim` is the dimension of the matrix.
"""
function companion(::Type{T}, v::AbstractVector) where {T}
    n = length(v)
    A = zeros(T, n, n)
    A[:, end] = v
    for i = 1:n-1
        A[i+1, i] = one(T)
    end
    A
end
companion(::Type{T}, n::Integer) where {T} = companion(T, [1:n;])
companion(args...) = companion(Float64, args...)
companion(::Type, args...) = throw(MethodError(companion, Tuple(args)))
