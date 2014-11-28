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
        lgn = iround(log2(n))
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


matrixdict = ["hilb" => hilb, "hadamard" => hadamard, 
              "cauchy" => cauchy, "circul" => circul,
              "dingdong" => dingdong, "frank" => frank,
              "invhilb" => invhilb, "forsythe" => forsythe];

matrixinfo = ["hilb" => "Hilbert matrix: 
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
              ];

matrixclass = ["symmetric" => ["hilb", "cauchy", "circul", "dingdong", 
                               "invhilb"],
               "inverse" => ["hilb", "hadamard", "cauchy", "invhilb", 
                             "forsythe"],
               "ill-cond" => ["hilb", "cauchy", "frank", "invhilb", 
                              "forsythe"],
               "pos-def" => ["hilb", "cauchy", "circul", "invhilb"],
               "eigen" =>   ["hadamard", "circul", "dingdong", "frank",
                             "eigen"],
               ];
