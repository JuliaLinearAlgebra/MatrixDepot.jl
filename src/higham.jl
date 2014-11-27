# Higham Test Matrices
export hilb, hadamard, cauchy, circul, dingdong, matrixdict, 
matrixclass, matrixinfo


#
# Hilbert Matrix
#
function hilb{T}(::Type{T}, m::Int, n::Int)
    # compute the Hilbert matrix
    H = zeros(T, m,n)
    for i =1:m, j = 1:n
        H[i,j] = one(T)/ (i + j -1)
    end
    return H
end
hilb{T}(::Type{T}, n::Int) = hilb(T, n, n)

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
    [D[i,j] = one(T)/2*(n-i-j+1.5) for i = 1:n, j= 1:n]
    return D
end

matrixdict = ["hilb" => hilb, "hadamard" => hadamard, 
              "cauchy" => cauchy, "circul" => circul,
              "dingdong" => dingdong];

matrixinfo = ["hilb" => "Hilbert matrix: 
              \n Input options: 
              \n (type), dim: the dimension of the matrix
              \n (type), row_dim, col_dim: the row and column dimension 
              \n ['inverse', 'ill-cond', 'symmetric', 'pos-def']",
              "hadamard" => "Hadamard matrix: 
              \n Input options: 
              \n (type), n::Int: the dimension of the matrix, n is a power of 2 
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
              \n (type), n:: Int: the dimension of the matrix.
              \n ['symmetric', 'eigen']",
              ];

matrixclass = ["symmetric" => ["hilb", "cauchy", "circul", "dingdong"],
               "inverse" => ["hilb", "hadamard", "cauchy"],
               "ill-cond" => ["hilb", "cauchy"],
               "pos-def" => ["hilb", "cauchy", "circul"],
               "eigen" =>   ["hadamard", "circul", "dingdong"],
               ];
