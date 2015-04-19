# Higham Test Matrices


#
# Hilbert Matrix
#
function hilb{T}(::Type{T}, m::Int, n::Int)
    # compute the Hilbert matrix
    H = zeros(T, m, n)
    for i =1:m, j = 1:n
        H[i,j] = one(T)/ (i + j -1)
    end
    return H
end
hilb{T}(::Type{T}, n::Int) = hilb(T, n, n)

#
# Inverse of Hilbert Matrix
#
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


#
# Hadamard Matrix
#
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


#
# Cauchy Matrix
#
function cauchy{T}(x::Vector{T}, y::Vector{T})    
    # Compute the cauchy matrix
    m = size(x,1)
    n = size(y,1)
    C = zeros(T, m, n)
    [C[i,j]= one(T)/(x[i] + y[j]) for i = 1:m, j = 1:n]
    return C
end

cauchy{T}(x::Vector{T}) = cauchy(x, x)
cauchy{T}(::Type{T}, k::Int) = cauchy(T[1:k;])

#
# Circul Matrix
#
function circul{T}(v::Vector{T}, w::Vector{T})
    # Compute the circul matrix
    # v: the first row of the matrix
    # w: the length of the vector is the dimension of the column
    m = length(v)
    n = length(w)
    C = zeros(T, m, n)
    [C[i,j] = v[1 + mod(j - i, n)] for i = 1:m, j = 1:n]
    return C
end

circul{T}(v::Vector{T}, n::Int) = circul(v, T[1:n])
circul{T}(v::Vector{T}) = circul(v, v)
circul{T}(::Type{T}, k::Int) = circul(T[1:k])

#
# Dingdong Matrix
#
function dingdong{T}(::Type{T}, n::Int)
    # Compute the dingdong matrix
    # type: data type
    # n: the dimension of the matrix
    D = zeros(T, n, n)
    [D[i,j] = one(T)/(2*(n-i-j+1.5)) for i = 1:n, j= 1:n]
    return D
end

# 
# Frank Matrix
#
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

#
# Forsythe Matrix
#
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
   
#
# Magic Matrix
#
function magic{T}(::Type{T}, n::Int)
    # Compute a magic square of order n
    # Learnt from Cleve Moler, Experiments with MATLAB, 2011
    if mod(n, 2) == 1
        # n is odd
        M = oddmagic(T, n)
    elseif mod(n, 4) == 0
        # n is doubly even
        if VERSION < v"0.4.0-dev+1827"
            a = ifloor(mod([1:n;], 4)/2)
        else
            a = floor(Integer, mod([1:n;], 4)/2)
        end
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

#
# Grcar Matrix
#
function grcar{T}(::Type{T}, n::Int, k::Int = 3)
    # Compute grcar matrix
    G = tril(triu(ones(T, n,n)), k) - diagm(ones(T, n-1), -1)
    return G
end

#
# Triw Matrix
#
function triw{T}(::Type{T}, m::Int, n::Int, alpha, k::Int)
    alpha = convert(T, alpha)
    A = tril(eye(T, m,n) + alpha * triu(ones(T, m,n), 1), k)
    return A
end
triw{T}(::Type{T}, n::Int) = triw(T, n, n,  -1, n-1)

#
# Moler Matrix
#
function moler{T}(::Type{T}, n::Int, alpha = -1.)
    alpha = convert(T, alpha)
    M = triw(T, n, n, alpha, n - 1)' * triw(T, n, n, alpha, n -1 ) 
    return M
end

#
# Pascal Matrix
#
function pascal{T}(::Type{T}, n::Int)
    P = zeros(T, n,n)
    [P[i,j] = binomial(i+j-2, j-1) for i = 1:n, j= 1:n]
    return P
end

# 
# Kahan Matrix
#
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

#
# Pei Matrix
#
function pei{T}(::Type{T}, n::Int, alpha = 1)
    alpha = convert(T, alpha)
    return alpha*eye(T, n, n) + ones(T, n, n)
end

# 
# Vandermonde Matrix
#
function vand{T}(p::Vector{T}, n::Int)
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
vand{T}(::Type{T}, n::Int) = vand(T[1:n;], n)
vand{T}(p::Vector{T}) = vand(p, length(p))

#
# Involutory Matrix
#
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

#
# Chebyshev spectral differetiation matrix.
#
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
        x[i] = cos (pi * (i - 1) / (n - 1))
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

#
# Lotkin Matrix
#
function lotkin{T}(::Type{T}, n::Int)
    A = hilb(T, n)
    A[1,:] = ones(T,n)'
    return A
end

#
# Clement Matrix
#
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

#
# Fiedler Matrix
#
function fiedler{T}(v::Vector{T})
    n = length(v)
    v = v[:].'
    A = ones(T, n) * v
    A = abs(A - A.') # nonconjugate transpose
end
fiedler{T}(::Type{T}, n::Int) = fiedler(T[1:n;])

#
# MIN[I,J] Matrix
#
function minij{T}(::Type{T}, n::Int)
    A = zeros(T, n, n)
    [A[i,j] = min(i,j) for i = 1:n, j = 1:n]
    return A
end

#
# Binomial Matrix
#
function binomialm{T}(::Type{T}, n::Int)
    # Mulitiple of involutory matrix
    L = Array(T, n, n)
    D = Diagonal((-2).^[0:n-1;])
    [L[i,j] = binomial(i-1, j-1) for i = 1:n, j = 1:n]
    U = L[n:-1:1, n:-1:1]
    return L*D*U
end

# 
# Tridiagonal Matrix
#
function tridiag{T}(x::Vector{T}, y::Vector{T}, z::Vector{T})
    return Tridiagonal(x,y,z)
end
# Toeplitz tridiagonal matrix
tridiag{T}(::Type{T}, n::Int, x::Int, y::Int, z::Int) = 
n == 1 ? y*ones(T,1,1) : 
         tridiag(x*ones(T, n-1), y*ones(T, n), z*ones(T, n-1))
tridiag{T}(::Type{T}, n::Int) = tridiag(T, n, -1, 2, -1) 
    
#
# Lehmer Matrix
#
function lehmer{T}(::Type{T}, n::Int)
    A = Array(T, n, n)
    [A[i,j] = min(i,j) / max(i,j) for i = 1:n, j = 1:n]
    return A
end

#
# Parter Matrix
#
function parter{T}(::Type{T}, n::Int)
    A = Array(T, n, n)
    [A[i,j] = one(T) / (i - j + 0.5) for i = 1:n, j = 1:n]
    return A
end

#
# Chow Matrix
#
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

#
# Random Correlation Matrix
#
# Reference: Numerically Stable Generation of Correlation Matrices
# and Their Factors. Philip Davies and Nicholas Higham, 
# BIT Numerical Mathematics, 2000, Vol 40. Issue 4, pp 640-651
#
# limitation: only generate matrix of type Float64
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

#
# Poisson Matrix
#
function poisson{T}(::Type{T}, n::Int)
    S = full(tridiag(T, n))
    A = speye(T, n)
    return kron(A,S) + kron(S,A)
end

#
# Neumann Matrix
#     
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

#
# Rosser Matrix
#
# References:
# 1. Rosser, Lanczos, Hestenes and Karush, J. Res. Natl. Bur. Stand. Vol. 47 (1951), pp291-297. 
# 2. Sylvester, Phil. Mag. ser.4 v.33-34 (1867), pp461
# For n = 8, a = 2, b = 1, the generated matrix is the test matrix used in reference 1.
#
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


#
# Matrix with application in sampling theory
#
function sampling{T}(x::Vector{T})
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
    return sampling(p)
end

#
# Wilkinson Matrix
#
function wilkinson{T}(::Type{T}, n::Int)
    if n == 1 # handle 1-d case
        return ones(T, 1, 1)
    end
    m = (n-1)/2
    A = Tridiagonal(ones(T,n-1), abs(T[-m:m;]), ones(T, n-1))
    return A
end

#
# Random matrix with element -1, 0, 1
#
function rando{T}(::Type{T}, m::Int, n::Int, k::Int)
    A = Array(T, m, n)    
    if k == 1
        copy!(A, floor(rand(m,n) + .5))
    elseif k == 2
        copy!(A, 2 * floor(rand(m,n) + .5) - one(T))
    elseif k == 3
        copy!(A, round(3 * rand(m,n) - 1.5))
    else
        error("invalid k value.")
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

#
# Random matrix with pre-assigned singular values
#
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
        error("invalid mode value.")
    end
    
    A = zeros(T, m, n)
    A[1:p, 1:p] = diagm(sigma)
    A = qmult!(A')
    A = qmult!(A')

    return A
end
randsvd{T}(::Type{T}, n::Int, kappa, mode) = randsvd(T, n, n, kappa, mode)
randsvd{T}(::Type{T}, n::Int, kappa) = randsvd(T, n, kappa, 3)
randsvd{T}(::Type{T}, n::Int) = randsvd(T, n, sqrt(1/eps()))

#
# Random orthogonal upper Hessenberg matrix
#
function rohess{T}(::Type{T}, n::Int)
    x = rand(n-1)*2*pi
    H = eye(T, n)
    H[n,n] = sign(randn())
    for i = n:-1:2
        theta = x[i-1]
        c = convert(T, cos(theta))
        s = convert(T, sin(theta))
        H[[i-1; i], :] = [c*H[i-1, :] + s*H[i, :]; -s*H[i-1, :] + c*H[i, :]]
    end
    return H
end

#
# KMS (Kac-Murdock-Szego) Toeplitz matrix
#

function kms{T}(::Type{T}, n::Int, rho)    
    A = Array(T, n, n)
    [A[i,j] = rho^(abs(i-j)) for i = 1:n, j = 1:n]
    if typeof(rho) <: Complex
        A = conj(tril(A, -1)) + triu(A)
    end
    return A
end
kms{T}(::Type{T}, n::Int) = kms(T, n, convert(T, 0.5))

#
# Wathen matrix
# A finite element matrix (sparse random entries).
# Reference:
# A.J. Wathen, Realistic eigenvalue bounds for the Galerkin
# mass matrix, IMA J, Numer. Anal., 7 (1987), pp. 449-457
#
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
                          "kms" => kms, "wathen" => wathen
                          );

matrixinfo = 
@compat Dict("hilb" => "Hilbert matrix: 
             \n Input options: 
             \n [type,] dim: the dimension of the matrix
             \n [type,] row_dim, col_dim: the row and column dimension 
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']",

             "invhilb" => "Inverse of Hilbert matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix
             \n ['inverse', 'ill-cond', 'symmetric','pos-def']",

             "hadamard" => "Hadamard matrix: 
             \n Input options: 
             \n [type,] dim: the dimension of the matrix, n is a power of 2 
             \n ['inverse', 'orthogonal', 'eigen']",

             "cauchy" => "Cauchy matrix: 
             \n Input options: 
             \n vec1, vec2: two vectors 
             \n vec: a vector
             \n [type,] dim: the dimension of the matrix
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']",

             "circul" => "Circul matrix: 
             \n Input options: 
             \n vec, n: a vector and the column dimension 
             \n vec: a vector 
             \n [type,] dim: the dimension of the matrix
             \n ['symmetric', 'pos-def', 'eigen']",

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
             \n ['ill-cond', 'eigen']",

             "forsythe" => "Forsythe matrix:
             \n Input options:
             \n [type,] n, alpha, lambda: n is the dimension of the matrix.
             alpha and lambda are scalars.
             \n [type,] n: alpha = sqrt(eps(type)) and lambda = 0.
             \n ['inverse', 'ill-cond', 'eigen']",

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
             \n ['inverse', 'ill-cond']",

             "moler" => "Moler Matrix:
             \n Input options:
             \n [type,] dim, alpha: dim is the dimension of the matrix,
             alpha is a scalar.
             \n [type,] dim: alpha = -1.
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']",

             "pascal" => "Pascal Matrix:
             \n Input options:
             \n [type,] dim: the dimension of the matrix.
             \n ['Inverse', 'ill-cond', 'symmetric', 'pos-def', 'eigen']",

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
             \n ['inverse', 'ill-cond', 'eigen']",

             "chebspec" => "Chebyshev spectral differentiation matrix:
             \n Input options:
             \n [type,] dim, k: dim is the dimension of the matrix and 
             k = 0 or 1.
             \n [type,] dim
             \n ['eigen']",

             "lotkin" => "Lotkin Matrix:
             \n Input options:
             \n [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'ill-cond', 'eigen']",

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
             \n ['inverse', 'symmetric', 'eigen']",

             "minij" => "MIN[I,J] Matrix:
             \n Input options:
             \n [type,] dim: dim is the dimension of the matrix.
             \n ['inverse', 'symmetric', 'pos-def', 'eigen']",

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
             \n ['inverse', 'symmetric', 'pos-def', 'eigen', 'sparse']",

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
             \n ['symmetric', 'eigen']",

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
             \n ['random']",

             "kms" => "Kac-Murdock-Szego Toeplitz Matrix:
             \n [type,] n, rho: n is the dimension of the matrix, rho is a 
             scalar such that A[i,j] = rho^(abs(i-j)).
             \n [type,] n: rho = 0.5
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']",

             "wathen" => "Wathen Matrix:
             \n [type,] nx, ny: the dimension of the matrix is equal to
             3 * nx * ny + 2 * nx * ny + 1.
             \n [type,] n: nx = ny = n. 
             \n ['symmetric', 'pos-def', 'eigen', 'random', 'sparse']
             \n Reference: A.J. Wathen, Realistic eigenvalue bounds for 
             the Galerkin mass matrix, IMA J. Numer. Anal., 7 (1987), 
             pp. 449-457."

             );

matrixclass = 
@compat Dict("symmetric" => ["hilb", "cauchy", "circul", "dingdong", 
                             "invhilb", "moler", "pascal", "pei", 
                             "clement", "fiedler", "minij", "tridiag",
                             "lehmer", "randcorr", "poisson", "wilkinson",
                              "kms", "wathen"],
             "inverse" => ["hilb", "hadamard", "cauchy", "invhilb", 
                           "forsythe", "magic", "triw", "moler", "pascal",
                           "kahan", "pei", "vand", "invol", "lotkin",
                           "clement", "fiedler", "minij", "tridiag", 
                           "lehmer", "poisson", "kms" ],
             "ill-cond" => ["hilb", "cauchy", "frank", "invhilb", 
                            "forsythe", "triw", "moler", "pascal",
                            "kahan","pei", "vand", "invol", "lotkin",
                            "tridiag", "rosser", "randsvd", "kms"],
             "pos-def" => ["hilb", "cauchy", "circul", "invhilb", 
                           "moler", "pascal", "pei", "minij", "tridiag",
                           "lehmer", "poisson", "kms", "wathen"],
             "eigen" =>   ["hadamard", "circul", "dingdong", "frank",
                           "forsythe", "grcar", "pascal", "invol","chebspec",
                           "lotkin", "clement", "fiedler", "minij",
                           "tridiag", "parter", "chow", "poisson", "neumann",
                           "rosser", "sampling", "wilkinson","wathen"],
             # minor properties
             "sparse" => ["poisson", "neumann", "wathen"],
             "random" => ["rosser", "rando", "randcorr", "randsvd", "rohess", 
                          "wathen"]
               );
