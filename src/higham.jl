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
    2^lgn != n && throw(ArgumentError("n must be positive integer and a power of 2"))
    
    H = reshape(T[1], 1, 1)
    for i = 1:lgn
        H = [[H H], [H -H]]
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
cauchy{T}(::Type{T}, k::Int) = cauchy(T[1:k])

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
    [C[i,j] = v[one(T) + mod(j - i, n)] for i = 1:m, j = 1:n]
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
    p = T[n:-1:1];
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
        a = ifloor(mod([1:n], 4)/2)
        B = broadcast(==, a', a)
        M = broadcast(+, T[1:n:n^2]',T[0:n-1])
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
        i = [1:p]
        k = div(n-2, 4)
        j = [[1:k], [(n-k+2) : n]]
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
    V = ones(T, m, n)
    for i = 1:n
        V[:,i] = p.^(i-1)
    end
    return V
end
vand{T}(::Type{T}, n::Int) = vand(T[1:n], n)
vand{T}(p::Vector{T}, v::Vector{T}) = vand(p, length(v))

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

matrixdict = @compat Dict("hilb" => hilb, "hadamard" => hadamard, 
                          "cauchy" => cauchy, "circul" => circul,
                          "dingdong" => dingdong, "frank" => frank,
                          "invhilb" => invhilb, "forsythe" => forsythe,
                          "magic" => magic, "grcar" => grcar,
                          "triw" => triw, "moler" => moler,
                          "pascal" => pascal, "kahan" => kahan,
                          "pei" => pei, "vand" => vand,
                          "invol" => invol,);

matrixinfo = 
@compat Dict("hilb" => "Hilbert matrix: 
             \n Input options: 
             \n (type), dim: the dimension of the matrix
             \n (type), row_dim, col_dim: the row and column dimension 
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']",
             "invhilb" => "Inverse of Hilbert matrix:
             \n Input options:
             \n (type), dim: the dimension of the matrix
             \n ['inverse', 'ill-cond', 'symmetric','pos-def']",
             "hadamard" => "Hadamard matrix: 
             \n Input options: 
             \n (type), dim: the dimension of the matrix, n is a power of 2 
             \n ['inverse', 'orthogonal', 'eigen']",
             "cauchy" => "Cauchy matrix: 
             \n Input options: 
             \n (type), vec1, vec2: two vectors 
             \n (type), vec: a vector
             \n (type), dim: the dimension of the matrix
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']",
             "circul" => "Circul matrix: 
             \n Input options: 
             \n (type), vec, col_dim: a vector and the column dimension 
             \n (type), vec: a vector 
             \n (type), dim: the dimension of the matrix
             \n ['symmetric', 'pos-def', 'eigen']",
             "dingdong" => "Dingdong matrix:
             \n Input options:
             \n (type), n: the dimension of the matrix.
             \n ['symmetric', 'eigen']",
             "frank" => "Frank matrix:
             \n Input options:
             \n (type), n, k: n is the dimension of the matrix, k = 0 or 1.
             If k = 1 the matrix reflect about the anti-diagonal.
             \n (type), n: n is the dimension of the matrix.
             \n ['ill-cond', 'eigen']",
             "forsythe" => "Forsythe matrix:
             \n Input options:
             \n (type), n, alpha, lambda: n is the dimension of the matrix.
             alpha and lambda are scalars.
             \n (type), n: alpha = sqrt(eps(type)) and lambda = 0.
             \n ['inverse', 'ill-cond', 'eigen']",
             "magic" => "Magic square matrix:
             \n Input options:
             \n (type), dim: the dimension of the matrix.
             \n ['inverse']",
             "grcar" => "Grcar Matrix:
             \n Input options:
             \n (type), dim, k: dim is the dimension of the matrix and 
             k is the number of superdiagonals.
             \n (type), dim: the dimension of the matrix.
             \n ['eigen']",
             "triw" => "Triw Matrix:
             \n Input options:
             \n (type), row_dim, col_dim, alpha, k: row_dim and col_dim 
             are row and column dimension of the matrix. alpha is a 
             scalar representing the entries on the superdiagonals. 
             k is the number superdiagonals.
             \n (type), dim 
             \n ['inverse', 'ill-cond']",
             "moler" => "Moler Matrix:
             \n Input options:
             \n (type), dim, alpha: dim is the dimension of the matrix,
             alpha is a scalar.
             \n (type), dim: alpha = -1.
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']",
             "pascal" => "Pascal Matrix:
             \n Input options:
             \n (type), dim: the dimension of the matrix.
             \n ['Inverse', 'ill-cond', 'symmetric', 'pos-def', 'eigen']",
             "kahan" => "Kahan Matrix:
             \n Input options:
             \n (type), m, n, theta, pert: m, n are the row and column 
             dimensions of the matrix. theta and pert are scalars.
             \n (type), dim, theta, pert: dim is the dimension of the matrix.
             \n (type), dim
             \n ['inverse', 'ill-cond']" ,
             "pei" => "Pei Matrix:
             \n Input options:
             \n (type), dim, alpha: dim is the dimension of the matrix.
             alpha is a scalar.
             \n (type), dim
             \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']",
             "vand" => "Vandermonde Matrix:
             \n Input options:
             \n vec, dim: vec is a vector, dim is the number of columns.
             \n vec
             \n (type), dim
             \n ['inverse', 'ill-cond']",
             "invol" => "Involutory Matrix:
             \n Input options:
             \n (type), dim: dim is the dimension of the matrix.
             \n ['inverse', 'ill-cond', 'eigen']",
             );

matrixclass = 
@compat Dict("symmetric" => ["hilb", "cauchy", "circul", "dingdong", 
                             "invhilb", "moler", "pascal", "pei", ],
             "inverse" => ["hilb", "hadamard", "cauchy", "invhilb", 
                           "forsythe", "magic", "triw", "moler", "pascal",
                           "kahan", "pei", "vand", "invol",],
             "ill-cond" => ["hilb", "cauchy", "frank", "invhilb", 
                            "forsythe", "triw", "moler", "pascal",
                            "kahan","pei", "vand", "invol",],
             "pos-def" => ["hilb", "cauchy", "circul", "invhilb", 
                           "moler", "pascal", "pei",],
             "eigen" =>   ["hadamard", "circul", "dingdong", "frank",
                           "forsythe", "grcar", "pascal", "invol",],
               );
