
"""
    mmread(filename|io)

Read `Matrixmarket` format file (extension `.mtx`) and return sparse or dense matrix.
Symmetric and Hermitian matrices use the corresponding wrapper types.
Patterns result in sparse matrics with element type `Bool`.
They may be converted to numerical types by multiplying with a number.
"""
function mmread(filename::AbstractString)
    open(filename, "r") do file
        mmread(file)
    end
end

using Mmap

const COORD = "coordinate"
const ARRAY = "array"
const MATRIX = "matrix"
const MATRIXM = "%%matrixmarket"
const COMPLEX = "complex"
const REAL = "real"
const INTEGER = "integer"
const PATTERN = "pattern"
const GENERAL = "general"
const SYMMETRIC = "symmetric"
const HERMITIAN = "hermitian"
const SKEW_SYMMETRIC = "skew-symmetric"

function mmread(file::IO)
    line = lowercase(readline(file))
    tokens = split(line)
    if tokens[1] != MATRIXM
        parserr(string("Matrixmarket: invalid header:", line))
    end
    line = readline(file)
    while length(strip(line)) == 0 || line[1] == '%'
        line = readline(file)
    end
    if tokens[2] == MATRIX
        mmread_matrix(file, line, tokens[3:end]...)
    else
        parserr(string("Matrixmarket: unsupported type: ", line))
    end
end

# mmap for regular files - else read
function getbytes(io::IOStream)
    isfile(io) ? Mmap.mmap(io, grow=false, shared=false) : read(io)
end
getbytes(io::IO) = read(io)

function mmread_matrix(file::IO, line, form, field, symm)
    FMAP = Dict(REAL => (3, Float64),
                COMPLEX => (4, ComplexF64),
                INTEGER => (3, Int64),
                PATTERN => (2, Bool))

    ty, T = get(FMAP, field) do
            parserr("Matrixmarket: unsupported field $field (only real/complex/pattern)")
    end

    SMAP = Dict(GENERAL => (1, 0, Any),
                SYMMETRIC => (0, 1, Symmetric),
                SKEW_SYMMETRIC => (1, 1, Array),
                HERMITIAN => (0, 1, Hermitian))

    p1, pc, wrapper = get(SMAP, symm) do
        parserr("Matrixmarket: unsupported symmetry $symm (general/symmetric/hermitian,skew-symmetric)")
    end

    if form == COORD
        m, n, nz = parseint(line)
        b = getbytes(file)
        rv = Vector{Int}(undef, nz)
        cv = Vector{Int}(undef, nz)
        vv = Vector{T}(undef, nz)
        parseloop!(Val(ty), b, rv, cv, vv)
        result = mksparse!(m, n, rv, cv, vv)
    elseif form == ARRAY
        m, n = parseint(line)
        b = getbytes(file)
        p = 1
        result = zeros(T, m, n)
        for c = 1:n
            for r = (pc*c+p1):m
                p, v = parsenext(T, b, p)
                result[r,c] = v
            end
        end
    else
        parserr("Matrixmarket: unsupported format '$form'")
    end

    wrap(result, wrapper)
end

wrap(result, ::Type{T}) where T<:Union{Symmetric,Hermitian} = T(result, :L)
wrap(result, ::Type{Array}) = result - mtranspose(result)
wrap(result, ::Type{Any}) = result

function parseloop!(::Val{4}, c::Vector{UInt8}, rv, cv, vv::Vector{T}) where T<:Complex
    nz = length(rv)
    R = real(T)
    p = 1
    for i = 1:nz
        p, rv[i] = parsenext(Int, c, p)
        p, cv[i] = parsenext(Int, c, p)
        p, r = parsenext(R, c, p)
        p, s = parsenext(R, c, p)
        vv[i] = r + s*im
    end
end

function parseloop!(::Val{3}, c::Vector{UInt8}, rv, cv, vv::Vector{T}) where T <:Real
    nz = length(rv)
    p = 1
    for i = 1:nz
        p, rv[i] = parsenext(Int, c, p)
        p, cv[i] = parsenext(Int, c, p)
        p, vv[i] = parsenext(T, c, p)
    end
end

function parseloop!(::Val{2}, c::Vector{UInt8}, rv, cv, vv::Vector{T}) where T<:Number

    nz = length(rv)
    p = 1
    for i = 1:nz
        p, rv[i] = parsenext(Int, c, p)
        p, cv[i] = parsenext(Int, c, p)
    end
    fill!(vv, T(1))
end

function parseint(line::AbstractString)
    tokens = split(line)
    parse.(Int, tokens)
end

"""
    mksparse(m, n, rowval, colval, nzval)
    mksparse!(m, n, rowval, colval, nzval)

Construct a `SparseMatrixCSC` of dimensions `(m,n)` from the data given in the
three input vectors of equal lengths.

`mksparse!` destoys the content of `colval`.

`A[rowval[i],colval[i]] == nzval[i] for i ∈ 1:length(nzval)`. All other entries are zero.
"""
mksparse(m, n, rv, cv, vv) = mksparse!(m, n, rv, copy(cv), vv)

function mksparse!(m::Integer, n::Integer, rv::AbstractVector{Ti}, cv::AbstractVector{Ti},
                  vv::AbstractVector{Tv}) where {Ti<:Integer,Tv<:Number}

    nz = length(rv)
    length(cv) == nz == length(vv) || argerr("all vectors need same length")
    micv, mcv = extrema(cv)
    mirv, mrv = extrema(rv)
    micv > 0 && mcv <= n || daterr("all column indices must be >= 1 and <= $n")
    mirv > 0 && mrv <= m || daterr("all row indices must be >= 1 and <= $m")
    sizeof(Ti) <= 8 || argerr("Index type greater 64 bits not supported")
    sr = count_ones(Ti(nextpow(2, mrv + 1) - 1))
    sh = count_zeros(Ti(nextpow(2, mcv + 1) - 1))
    sr < sh || argerr("combined size of column and row indices exceeds $Ti size")
    # Ti must be able to keep (nz + 1)
    (nz + 1) % Ti != nz + 1 && argerr("Ti($Ti) cannot store nz($nz)")
    # compress row, col into Ti
    # t = UInt64[]
    # push!(t, time_ns())
    colptr = zeros(Ti, n+1)
    colptr[1] = 1
    # push!(t, time_ns())
    @inbounds for i = 1:nz
        cvi = cv[i]
        cv[i] = cvi << sr | rv[i]
        colptr[cvi+1] += 1
    end
    # push!(t, time_ns())
    cumsum!(colptr, colptr)
    # push!(t, time_ns())
    p = specialsort(cv, sr)
    # push!(t, time_ns())
    # println("times: $(diff(t) ./ 1e6) ms")
    SparseMatrixCSC{Tv,Ti}(m, n, colptr, rv[p], vv[p])
end

function specialsort(cv::Vector{Int}, sr::Int)
    nz = length(cv)
    if nz > 10000
        x = isqrt(nz) << sr
        p = sortperm(cv .÷ x)
        sortperm!(p, cv, initialized=true)
    else
        sortperm(cv)
    end
end

"""
    colval(A)
reconstruct colum-indices from colptr of `SparseMatrixCSC`.
"""
function colval(A::SparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    nz = nnz(A)
    cv = Vector{Ti}(undef, nz)
    colptr = A.colptr
    for j = 1:A.n
        for i = colptr[j]:colptr[j+1]-1
            cv[i] = j
        end
    end
    cv
end

"""
    mtranspose(A)
Materialized transpose of a matrix
"""
mtranspose(A::SparseMatrixCSC) = mksparse!(A.n, A.m, colval(A), copy(A.rowval), copy(A.nzval))
mtranspose(A::Matrix) = Matrix(transpose(A))
mtranspose(A) = transpose(A)

"""
    madjoint(A)
Materialized adjoint of sparse Matrix
"""
madjoint(A) = conj!(mtranspose(A))

"""
    mmreadcomment(filename)
return info comment strings for MatrixMarket format files
"""
function mmreadcomment(filename::AbstractString)
    io = IOBuffer()
    open(filename,"r") do mmfile
        line = readline(mmfile)
        while startswith(line, '%') || isempty(strip(line))
            println(io, line)
            line = readline(mmfile)
        end
        println(io, line)
    end
    String(take!(io))
end

"""
    mmreadheader(filename)
Read header information from mtx file.
"""
function mmreadheader(file::AbstractString)
    if isfile(file)
        open(file) do io
            line = lowercase(readline(io))
            token = split(line)
            if length(token) >= 4 &&
                token[1] == MATRIXM &&
                token[2] == MATRIX &&
                token[3] in [COORD, ARRAY]

                hdr = Dict{Symbol,Any}()
                field = :none
                while startswith(line, '%') || isempty(strip(line))
                    field = push_hdr!(hdr, line, field) 
                    line = readline(io)
                end
                res = try parseint(line) catch; [] end
                if length(res) != (token[3] == COORD ? 3 : 2)
                    daterr("MatrixMarket file '$file' invalid sizes: '$line'")
                end
                hdr[:m] = res[1]
                hdr[:n] = res[2]
                length(res) >= 3 && (hdr[:nz] = res[3])
                hdr[:format] = token[3]
                hdr[:field] = token[4]
                hdr[:symmetry] = token[5]
                if haskey(hdr, :notes)
                    hdr[:notes] = join(wordlist(String(take!(hdr[:notes]))), ' ')
                end
                if haskey(hdr, :date)
                    val = hdr[:date]
                    hdr[:date] = isempty(val) ? 0 : parse(Int, hdr[:date])
                end
                hdr
            else
                daterr("file '$file' is not a MatrixMarket file")
            end
        end
    else
        nothing
    end
end

"""
    wordlist(string)

Separate words is string by spaces and delimiters.
Return list of unique words, which can be used as keywords.
"""
function wordlist(s::AbstractString)
    list = unique!(split(s, r"[][\s(){}`\"'*]", keepempty = false))
    list = replace.(list, Ref(r"[.:,;']$" => ""))
    # remove all lowercase words with less than ... chars
    list = filter!(x->!(length(x)<4 && all(islowercase.(collect(x))) || length(x) < 2), list)
    unique!(list)
end

function push_hdr!(hdr, line::AbstractString, field::Symbol)
    isempty(strip(line)) && return field
    if field == :notes
        if !startswith(line, "%---")
            println(get!(hdr, field) do; IOBuffer() end, strip(line[2:end]))
        end
        return field
    end
    reg = r"^% *([^:[]+): *(.*)$"
    regtitle = r"^% *\[([^]]*)]"
    if (m = match(reg, line)) != nothing
        s = Symbol(m[1])
        if s in (:name, :kind, :ed, :fields, :author, :date)
            field = s
            value = strip(m[2])
            hdr[field] = value
        elseif s == :notes
            field = s
        end
    elseif (m = match(regtitle, line)) != nothing
        field = :title
        value = m[1]
        hdr[field] = value
    end
    field
end

### parsing decimal integers and floats
function _parsenext(v::Vector{UInt8}, p1::Int)
    iaccu = Unsigned(0)
    daccu = Unsigned(0)
    eaccu = 0
    df = 0
    sig = 0
    esig = 0
    i = p1
    n = length(v)
    c = v[i]
    while c == 0x20 || c == 0x0a || c == 0x0d || c == 0x09
        c = v[i += 1]
    end
    n0 = i
    if c == UInt8('+')
        sig = 1
        c = v[i += 1]
    elseif c == UInt8('-')
        sig = -1
        c = v[i += 1]
    end
    ne = i
    while 0x30 <= c <= 0x39
        iaccu = iaccu * 10 + c - 0x30
        c = v[i += 1]
    end
    if c == UInt8('.')
        c = v[i += 1]
        d0 = i
        while 0x30 <= c <= 0x39
            daccu = daccu * 10 + c - 0x30
            c = v[i += 1]
        end
        df = i - d0
        ne += 1
    end
    i > ne || parserr("Invalid decimal number: '$(String(v[n0:min(end,n0+5)]))'")
    if i > ne && ( c == UInt8('e') || c == UInt8('E') )
        c = v[i += 1]
        if c == UInt8('+')
            esig = 1
            c = v[i += 1]
        elseif c == UInt8('-')
            esig = -1
            c = v[i += 1]
        end
        while 0x30 <= c <= 0x39
            eaccu = eaccu * 10 + c - 0x30
            c = v[i += 1]
        end
        if esig < 0
            eaccu = -eaccu
        end
    end
    i == n0 && parserr("No decimal number found: '$(String(v[n0:min(end,n0+5)]))'")
    i, iaccu, daccu, eaccu, sig, df
end

function parsenext(T::Type{<:Signed}, v::Vector{UInt8}, p1::Int)
    i, iaccu, daccu, eaccu, sig, df = _parsenext(v, p1)
    daccu == 0 && eaccu == 0 && df == 0 || error("1")
    i, T(iaccu) * ifelse(sig < 0, T(-1), T(1))
end
function parsenext(T::Type{<:Unsigned}, v::Vector{UInt8}, p1::Int)
    i, iaccu, daccu, eaccu, sig, df = _parsenext(v, p1)
    daccu == 0 && eaccu == 0 && sig >= 0 && df == 0 || error("2")
    i, T(iaccu)
end
function parsenext(T::Type{<:AbstractFloat}, v::Vector{UInt8}, p1::Int)
    i, iaccu, daccu, eaccu, sig, df = _parsenext(v, p1)
    f = exp10(eaccu) * (T(daccu) / T(exp10(df)) + T(iaccu))
    i, (sig < 0 ? -f : f)
end
function parsenext(T::Type{<:Complex}, c, p)
    R = real(T)
    p, r = parsenext(R, c, p)
    p, s = parsenext(R, c, p)
    p, r + s*im
end

