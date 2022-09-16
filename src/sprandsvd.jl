using SparseArrays
using DataStructures
export SparseMatrixD

struct SparseMatrixD{T,S} <: AbstractSparseMatrix{T,S}
    m::Int
    n::Int
    cols::Vector{SortedDict{Int,T}}
    function SparseMatrixD(A::AbstractSparseMatrix{T,S}) where {T,S}
        m, n = size(A)
        cols = Vector{SortedDict{S,T}}(undef, n)
        rows = rowvals(A)
        vals = nonzeros(A)
        @inbounds for j in axes(A, 2)
            cols[j] = colj = SortedDict{S,T}()
            for k in nzrange(A, j)
                v = vals[k]
                if !iszero(v)
                    colj[rows[k]] = v
                end
            end
        end
        new{T,S}(m, n, cols)
    end
end

Base.size(S::SparseMatrixD) = (S.m, S.n)
SparseArrays.nnz(S::SparseMatrixD) = sum(length.(S.cols))

function Base.getindex(S::SparseMatrixD{T}, i::Integer, j::Integer) where {T}
    @boundscheck checkbounds(S, i, j)
    colj = S.cols[j]
    get(colj, i, zero(T))
end

function Base.setindex!(S::SparseMatrixD{T}, v::Any, i::Integer, j::Integer) where {T}
    @boundscheck checkbounds(S, i, j)
    colj = S.cols[j]
    colj[i] = v
end

function Base.show(io::IO, m::MIME{Symbol("text/plain")}, S::SparseMatrixD)
    show(io, m, SparseMatrixCSC(S))
end

function SparseArrays.SparseMatrixCSC(A::SparseMatrixD{T,S}) where {T,S}
    m, n = size(A)
    nz = nnz(A)
    colptr = Vector{S}(undef, n + 1)
    rowval = Vector{S}(undef, nz)
    nzval = Vector{T}(undef, nz)
    k = 1
    colptr[1] = k
    @inbounds for j in axes(A, 2)
        colj = A.cols[j]
        for (i, v) in colj
            if !iszero(v)
                @assert k <= nz
                rowval[k] = i
                nzval[k] = v
                k += 1
            end
        end
        colptr[j+1] = k
    end
    if nz != k - 1
        resize!(rowval, k - 1)
        resize!(nzval, k - 1)
    end
    SparseMatrixCSC{T,S}(m, n, colptr, rowval, nzval)
end

"""
Random Sparse Matrix with Pre-assigned Singular Values
======================================================
*Input options:*

+ [type,] row_dim, col_dim, p, sigma::Vector
  `sigma` is the vector of singular values
  `p` is the target density of the matrix (0 < p < 1)

+ [type,] row_dim, col_dim, p, κ, mode: `row_dim` and `col_dim`
    are the row and column dimensions.
  `κ` is the condition number of the matrix.
  `mode = 1`: one large singular value.
  `mode = 2`: one small singular value.
  `mode = 3`: geometrically distributed singular values.
  `mode = 4`: arithmetrically distributed singular values.
  `mode = 5`: random singular values with  unif. dist. logarithm.

+ [type,] dim, p, κ, mode: `row_dim = col_dim = dim`.

+ [type,] dim, p, κ: `mode = 3`.

+ [type,] dim, p: `κ = sqrt(1/eps())`, `mode = 3`.

*Groups:* ["illcond", "random"]

*References:*

Similar to `sprandn`

https://de.mathworks.com/help/matlab/ref/sprand.html
"""
function sprandsvd(::Type{T}, m::Integer, n::Integer, p::AbstractFloat,
    κ::AbstractFloat=sqrt(1 / eps(T)), mode::Integer=3; rng=Random.default_rng()) where {T}

    sigma = make_sigma(rng, min(m, n), (real(T))(κ), mode)
    sprandsvd(T, m, n, p, sigma; rng)
end

function sprandsvd(::Type{T}, n::Integer, p::AbstractFloat,
    κ::AbstractFloat, mode::Integer=3; kw...) where {T}

    sprandsvd(T, n, n, p, κ, mode; kw...)
end
function sprandsvd(::Type{T}, n::Integer, p::AbstractFloat; kw...) where {T}
    sprandsvd(T, n, n, p; kw...)
end
sprandsvd(args...; kw...) = sprandsvd(Float64, args...; kw...)
sprandsvd(::Type, args...; kw...) = throw(MethodError(sprandsvd, Tuple(args)))

function sprandsvd(::Type{T}, m::Integer, n::Integer, p::AbstractFloat,
    sigma::AbstractVector; rng=Random.default_rng()) where {T}

    rng::AbstractRNG
    0 <= p <= 1 || throw(ArgumentError("$p not in [0,1]"))
    m >= 0 && n >= 0 || throw(ArgumentError("invalid Array dimensions"))
    v = sigma
    mn = min(m, n, length(v))
    R = spzeros(T, m, n)
    (mn == 0 || iszero(v)) && return R
    S = SparseMatrixD(R)
    pm = randperm(rng, m)
    pn = randperm(rng, n)
    rows = BitSet()
    cols = BitSet()
    for i in 1:mn
        S[pm[i], pn[i]] = v[i]
        push!(rows, pm[i])
        push!(cols, pn[i])
    end
    m <= 1 && n <= 1 && return SparseMatrixCSC(S)
    nnzs = nnz(S)
    nnzm = max(Int(floor(p * m * n)), nnzs)
    cl = 1
    cr = 1
    @inbounds while nnzs < nnzm
        rowx = (m - 1) * length(rows) * cr
        colx = (n - 1) * length(cols) * cl
        if rand(rng, 1:rowx+colx) <= rowx
            i = rand(rng, rows)
            j = rand(rng, 1:m-1)
            j = ifelse(j < i, j, j + 1)
            g11, g12, g21, g22 = randu(rng, T)
            leftmul!(S, g11, g12, g21, g22, i, j)
            push!(rows, j)
            cl += 1
        else
            i = rand(rng, cols)
            j = rand(rng, 1:n-1)
            j = ifelse(j < i, j, j + 1)
            g11, g12, g21, g22 = randu(rng, T)
            rightmul!(S, g11, g12, g21, g22, i, j)
            push!(cols, j)
            cr += 1
        end
        nnzs = nnz(S)
    end
    return SparseMatrixCSC(S)
end

function randu(rng::AbstractRNG, ::Type{T}) where {T<:AbstractFloat}
    s, c = sincospi(rand(rng, T) * 2)
    @inbounds t = (-1, 1)[rand(rng, 1:2)]
    c, s, -s * t, c * t
end
function randu(rng::AbstractRNG, ::Type{Complex{T}}) where {T<:AbstractFloat}
    r = rand(rng, T, 4)
    s, c = sincospi(r[1] * 0.5)
    a = exp.(r[2:4] .* pi * 2 * im)
    @inbounds (c * a[1], s * a[2], -s * a[1] * a[3], c * a[2] * a[3])
end

function leftmul!(S::SparseMatrixD, g11, g12, g21, g22, i, j)
    @boundscheck checkbounds(S, i, 1)
    @boundscheck checkbounds(S, j, 1)
    @inbounds for k = axes(S, 2)
        si = S[i, k]
        sj = S[j, k]
        iszero(si) && iszero(sj) && continue
        S[i, k] = g11 * si + g12 * sj
        S[j, k] = g21 * si + g22 * sj
    end
    S
end
function rightmul!(S::SparseMatrixD{T}, g11, g12, g21, g22, i, j) where {T}
    @boundscheck checkbounds(S, 1, i)
    @boundscheck checkbounds(S, 1, j)
    @inbounds ci = S.cols[i]
    @inbounds cj = S.cols[j]
    @inbounds for k = axes(S, 1)
        si = get(ci, k, zero(T))
        sj = get(cj, k, zero(T))
        iszero(si) && iszero(sj) && continue
        ci[k] = g11 * si + g21 * sj
        cj[k] = g12 * si + g22 * sj
    end
    S
end

"""
Random Symmetric Sparse Matrix with Pre-assigned Singular Values
================================================================
*Input options:*

+ [type,] dim, p, sigma::Vector
  `sigma` is the vector of (real) eigenvalues
  `p` is the target density of the matrix (0 < p < 1)

+ [type,] dim, p, κ, mode: `dim` is the size of the square matrix.
  `κ` is the condition number of the matrix.
  `mode = 1`: one large eigenvalue.
  `mode = 2`: one small eigenvalue.
  `mode = 3`: geometrically distributed singular values.
  `mode = 4`: arithmetrically distributed singular values.
  `mode = 5`: random singular values with  unif. dist. logarithm.

+ [type,] dim, p, κ, mode

+ [type,] dim, p, κ: `mode = 3`.

+ [type,] dim, p: `κ = sqrt(1/eps())`, `mode = 3`.

*Groups:* ["illcond", "random", "symmetric"]

*References:*

Similar to `sprandn`

https://de.mathworks.com/help/matlab/ref/sprand.html
"""
function sprandsym(::Type{T}, n::Integer, p::AbstractFloat,
    κ::AbstractFloat, mode::Integer=3; rng=Random.default_rng()) where {T}

    sigma = make_sigma(rng, n, (real(T))(κ), mode)
    sprandsym(T, n, p, sigma; rng)
end
sprandsym(args...; kw...) = sprandsym(Float64, args...; kw...)
sprandsym(::Type, args...; kw...) = throw(MethodError(sprandsym, Tuple(args)))

function sprandsym(::Type{T}, n::Integer, p::AbstractFloat,
    sigma::AbstractVector{<:Real}; rng=Random.default_rng()) where {T}

    rng::AbstractRNG
    0 <= p <= 1 || throw(ArgumentError("$p not in [0,1]"))
    n >= 0 || throw(ArgumentError("invalid Array dimensions"))
    v = sigma
    mn = min(n, length(v))
    R = spzeros(T, n, n)
    (mn == 0 || iszero(v)) && return R
    S = SparseMatrixD(R)
    pn = randperm(rng, n)
    rows = BitSet()
    for i in 1:mn
        S[pn[i], pn[i]] = v[i]
        push!(rows, pn[i])
    end
    nnzs = nnz(S)
    nnzm = max(Int(floor(p * n * n)), nnzs)
    @inbounds while nnzs < nnzm
        i = rand(rng, rows)
        j = rand(rng, 1:n-1)
        j = ifelse(j < i, j, j + 1)
        g11, g12, g21, g22 = randu(rng, T)
        leftmul!(S, g11, g12, g21, g22, i, j)
        rightmul!(S, conj(g11), conj(g21), conj(g12), conj(g22), i, j)
        push!(rows, j)
        nnzs = nnz(S)
    end
    return SparseMatrixCSC(symmetrize!(S))
end
