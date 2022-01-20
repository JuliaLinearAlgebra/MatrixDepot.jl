
# Exception

struct DataError <:Exception
    msg::String
end

## MatrixMarket header

abstract type MMProperty end
abstract type MMObject <: MMProperty end
struct MMObjectMatrix <: MMObject end

abstract type MMFormat <: MMProperty end
struct MMFormatCoordinate <: MMFormat end
struct MMFormatArray <: MMFormat end

abstract type MMField <:MMProperty end
struct MMFieldReal <: MMField end
struct MMFieldComplex <: MMField end
struct MMFieldInteger <: MMField end
struct MMFieldPattern <: MMField end

abstract type MMSymmetry <: MMProperty end
struct MMSymmetryGeneral <: MMSymmetry end
struct MMSymmetrySymmetric <: MMSymmetry end
struct MMSymmetrySkewSymmetric <: MMSymmetry end
struct MMSymmetryHermitian <: MMSymmetry end

struct MMProperties
    object::MMObject
    format::MMFormat
    field::MMField
    symmetry::MMSymmetry
end

mutable struct MetaInfo
    m::Int          # number of rows
    n::Int          # number of columns
    nnz::Int        # number of nonzeros in sparse matrix
    dnz::Int        # number of nonzero data items in file (half of nnz for symmetric)
    kind::String    # category of problem
    date::Int # the Julian year number e.g. 2018
    title::AbstractString   #
    author::AbstractString  #
    ed::AbstractString      # editor
    fields::AbstractString  # names of metadata and files like :A, :b, :x
    notes::AbstractString   # list of relevant words in notes
    # the following data are obtained from suite sparse index file ss_index.mat
    nnzdiag::Int    # number of zeros in diagonal
    pattern_symmetry::Float64   # [0..1] symmetry of non-zero element patterns
    numerical_symmetry::Float64 # [0..1] symmetry regarding A[i,j] == A[j,i]
    posdef::Bool    # Matrix is positive definite 1 (and symmetric/hermitian)
    isND::Bool      # problem comes form 2D/3D discretization
    isGraph::Bool   # problem is best considered as a graph than a system of equations
    cholcand::Bool  # matrix is symmetric (Hermitian if complex) and diagonal is positive
    amd_lnz::Int    # -2 if not square, -1 if not computed, else estimation for factorization
    amd_vnz::Int    # upper bound on the number of entries in L of LU factorization
    amd_rnz::Int    # upper bound on the number of entries in U in LU factorization
    amd_flops::Int  # floating operation count for Cholesky factorization of symm. pattern
    ncc::Int        # number of strongly connected components of graph/bipartite graph of A
    nblocks::Int    # number of blocks in dmperm (MATLAB)
    sprank::Int     # structural rank of A, max number of nonzero entries which can be
                    # permuted to diagonal (dmperm)
    # the following block of data is not well documented but delivered by suite sparse
    lowerbandwidth::Int
    upperbandwidth::Int
    rcm_lowerbandwidth::Int
    rcm_upperbandwidth::Int
    xmin::ComplexF64
    xmax::ComplexF64
    # the following data are derived from the singular value decomposition of A
    svdok::Bool      # singular value decomposition was successfull - data available
    norm::Float64   # maximal SV
    minsv::Float64  # minimal positive SV
    cond::Float64   # condition number in 2-norm (norm/minsv)
    rank::Int       # numerical rank (num. of SV > tol = max(m,n) * eps(norm)
    nullspace::Int  # dimension of null space of A = min(m,n) - rank
    svgap::Float64  # SV[rank]/SV[rank+1] if length(sv), else Inf
    svdstatus::String # "" if no SV data available, "OK" if converged, or a warning text
    svdhow::String  # "text describing how svd values was calculated"
    sv::Vector{Float64} # singular values > 0
end

## Matrix objects

struct RemoteParameters
    site::String
    dataurl::String
    indexurl::String
    indexgrep::String
    scan::Tuple
    extension::String
end
struct RemoteParametersNew
    site::String
    dataurl::String
    indexurl::String
    statsdb::String
    extension::String
end

abstract type RemoteType end

struct SSRemoteType <: RemoteType
    params::RemoteParametersNew
end

struct MMRemoteType <: RemoteType
    params::RemoteParameters
end

abstract type MatrixData end

struct RemoteMatrixData{T<:RemoteType} <:MatrixData
    name::AbstractString
    id::Int
    header::MetaInfo
    properties::Ref{Union{<:MMProperties,Nothing}}
    metadata::Vector{AbstractString}
    function RemoteMatrixData{T}(name, id::Integer, hdr::MetaInfo) where T
        properties = Ref{Union{<:MMProperties,Nothing}}(nothing)
        new(name, id, hdr, properties, AbstractString[])
    end
end

struct GeneratedMatrixData{N} <:MatrixData
    name::AbstractString
    id::Int
    func::Function
end

struct MatrixDatabase
    data::Dict{AbstractString,MatrixData}
    aliases::Dict{AbstractString,AbstractString}
    MatrixDatabase() = new(Dict{AbstractString,MatrixData}(),
                           Dict{AbstractString,AbstractString}())
    MatrixDatabase(data, aliases) = new(data, aliases)
end
Base.show(io::IO, db::MatrixDatabase) = print(io, "MatrixDatabase(", length(db.data), ")")

# Patterns

const IntOrVec = Union{Integer,AbstractVector{<:Integer}}
const AliasArgs = Union{Integer,Colon,AbstractVector{<:Union{Integer,AbstractRange{<:Integer}}}}

struct Alias{T,D}
    data::D
    Alias{T}(d::D) where {T<:MatrixData,D} = new{T,D}(d)
    function Alias{T}(d...) where {T<:MatrixData}
        v = getindex.(convertalias(collect(d))) # strip Ref wrap if necessary
        Alias{T}(besttype(v))
    end
end

convertalias(a::Integer) = a
convertalias(a::Colon) = Ref(a)
convertalias(a::AbstractRange{<:Integer}) = Ref(a)
convertalias(a::AbstractVector) = begin x = (vcat(convertalias.(a)...)); length(x) == 1 ? x[1] : x end

besttype(a::AbstractVector{T}) where T<:Integer = a
besttype(a::AbstractVector{T}) where T<:AbstractRange = a
besttype(a::AbstractVector{Any}) = (:) in a ? (:) : Vector{Union{Integer,AbstractRange}}(a)

struct Alternate{T<:RemoteType,P}
    pattern::P
end

abstract type AbstractNot end

const Pattern = Union{Function,AbstractString,Regex,Symbol,Alias,
                      Alternate,AbstractVector,Tuple,AbstractNot}

struct Not{T<:Pattern} <:AbstractNot
    pattern::T
end

"""
    MatrixDescriptor{T<:MatrixData}

Access object which is created by a call to `mdopen(::Pattern)`.
Several user functions allow to access data according to the unique pattern.
Keeps data cache for exactly the same calling arguments (in the case of
generated objects). For remote objects the data cache keeps all available
metadata.
"""
struct MatrixDescriptor{T<:MatrixData}
    data::T
    args::Tuple
    cache::Union{Ref,Dict}
end
function MatrixDescriptor(data::T) where T<:RemoteMatrixData
    MatrixDescriptor{T}(data, (), Dict())
end

function MatrixDescriptor(data::T, args...) where T<:GeneratedMatrixData
    MatrixDescriptor{T}(data, deepcopy(args), Ref{Any}(nothing))
end

struct Auxiliar{T<:MatrixDescriptor}
    md::T
end
# essential functions of the types

function Base.show(io::IO, data::RemoteMatrixData)
    hd = data.header
    prop = data.properties[]
    print(io, "(")
    print(io, prop === nothing ? "" : string(prop, " ") )
    print(io, "$(data.name)($(aliasname(data))) ")
    nnz = hd.nnz == hd.dnz ? "$(hd.nnz)" : "$(hd.nnz)/$(hd.dnz)"
    print(io, " $(hd.m)x$(hd.n)($nnz) ")
    print(io, data.date != 0 ? data.date : "")
    meta = join(metastring.(data.name, metadata(data)), ", ")
    n = length(meta)
    if n > 40
        meta = string(meta[1:17], " ... ", meta[end-17:end])
    end
    print(io, " [", meta, "]")
    print(io, " '", data.kind, "'")
    print(io, " [", data.title, "]")
    print(io, ")")
end

function Base.show(io::IO, mdesc::MatrixDescriptor)
    show(io, mdesc.data)
    show(io, mdesc.args)
end

function Base.push!(db::MatrixDatabase, data::MatrixData)
    key = data.name
    db.data[data.name] = data
    db.aliases[aliasname(data)] = key
end

"""
    aliasname(data::MatrixData)
    aliasname(Type{<:MatrixData, id::Integer)

return alias name derived from integer id
"""
aliasname(::Type{RemoteMatrixData{SSRemoteType}}, i::Integer) = string('#', i)
aliasname(::Type{RemoteMatrixData{MMRemoteType}}, i::Integer) = string('#', 'M', i)
aliasname(::Type{GeneratedMatrixData{N}}, i::Integer) where N = string('#', N, i)
function aliasname(T::Type{<:MatrixData}, r::AbstractVector)
    aliasname.(T, [ x for x in Iterators.flatten(r) if x > 0])
end
aliasname(T::Type{<:MatrixData}, ::Colon) = aliasname(T, 0)
aliasname(data::MatrixData) = aliasname(typeof(data), data.id)
aliasname(ali::Alias{T,D}) where {T,D} = aliasname(T, ali.data)

Base.show(io::IO, a::Alias) = print(io, aliasname(a))

import Base: ==
==(a::S, b::T) where {R,S<:Alias{R},T<:Alias{R}} = aliasname(a) == aliasname(b)

Base.empty!(db::MatrixDatabase) = (empty!(db.aliases); empty!(db.data))

localindex(::SSRemoteType) = abspath(data_dir(), "ss_index.mat")
localindex(::MMRemoteType) = abspath(data_dir(), "mm_matrices.html")

function sitename(s::AbstractString, t::AbstractString)
    p = split(s, '/')
    string((length(p) >= 3 ? p[3] : s), " with ", t, " index")
end

remote_name(s::SSRemoteType) = sitename(s.params.site, "MAT")
remote_name(s::MMRemoteType) = sitename(s.params.site, "HTML")

localbase(::Type{SSRemoteType}) = abspath(data_dir(), "uf")
localbase(::Type{MMRemoteType}) = abspath(data_dir(), "mm")

localdir(data::RemoteMatrixData{T}) where T = abspath(localbase(T), filename(data))
localdir(data::MatrixData) = nothing

indexurl(remote::RemoteType) = remote.params.indexurl
dataurl(remote::RemoteType) = remote.params.dataurl
dataurl(::Type{T}) where T<:RemoteType = dataurl(preferred(T))
extension(::Type{T}) where T<:RemoteType = extension(preferred(T))
extension(remote::RemoteType) = remote.params.extension
dataurl(data::RemoteMatrixData{T}) where T = join((dataurl(T), data.name), '/') * extension(T)
siteurl(data::RemoteMatrixData{T}) where T = preferred(T).params.site
filename(data::RemoteMatrixData) = joinpath(split(data.name, '/')...)
localfile(data::RemoteMatrixData{<:SSRemoteType}) = localdir(data) * ".tar.gz"
localfile(data::RemoteMatrixData{MMRemoteType}) = localdir(data) * ".mtx.gz"
function matrixfile(data::RemoteMatrixData{<:SSRemoteType})
    joinpath(localdir(data), rsplit(data.name, '/', limit=2)[end] * ".mtx")
end
matrixfile(data::RemoteMatrixData{MMRemoteType}) = localdir(data) * ".mtx"
function matrixinfofile(data::RemoteMatrixData)
    string(rsplit(matrixfile(data), '.', limit=2)[1], ".info")
end

## MatrixMarket header

mm_property_name(::MMObjectMatrix) = "matrix"
mm_property_name(::MMFormatCoordinate) = "coordinate"
mm_property_name(::MMFormatArray) = "array"
mm_property_name(::MMFieldReal) = "real"
mm_property_name(::MMFieldComplex) = "complex"
mm_property_name(::MMFieldInteger) = "integer"
mm_property_name(::MMFieldPattern) = "pattern"
mm_property_name(::MMSymmetryGeneral) = "general"
mm_property_name(::MMSymmetrySymmetric) = "symmetric"
mm_property_name(::MMSymmetrySkewSymmetric) = "skew-symmetric"
mm_property_name(::MMSymmetryHermitian) = "hermitian"

mm_property_type(::MMFieldReal) = Float64
mm_property_type(::MMFieldComplex) = ComplexF64
mm_property_type(::MMFieldInteger) = Integer
mm_property_type(::MMFieldPattern) = Bool

MM_NAME_TO_PROP = Dict{String,MMProperty}(mm_property_name(x) => x for x in (
    MMObjectMatrix(),
    MMFormatCoordinate(), MMFormatArray(),
    MMFieldReal(), MMFieldComplex(), MMFieldInteger(), MMFieldPattern(),
    MMSymmetryGeneral(), MMSymmetrySymmetric(), MMSymmetrySkewSymmetric(),
    MMSymmetryHermitian())
)

MMProperties() = MMProperties("matrix", "coordinate", "real", "general")
function MMProperties(args::AbstractString...)
    prop(x) = MM_NAME_TO_PROP[lowercase(x)]
    MMProperties(prop.(args)...)
end
Base.show(io::IO, pr::MMProperty) = print(io, mm_property_name(pr))
function Base.show(io::IO, mp::MMProperties)
    mms(pr::MMProperty) = mm_short_name(pr)
    print(io, #= mms(mp.object), mms(mp.format),=# mms(mp.field), mms(mp.symmetry))
end
mm_short_name(pr::MMSymmetrySkewSymmetric) = "K"
mm_short_name(pr::MMProperty) = uppercase(first(mm_property_name(pr)))

"""
    flattenPattern(p::Pattern)
return the vector of all elementary patterns, contained in the pattern.
"""
flatten_pattern(p::Pattern) = collect(Set(_flatten_pattern(p)))
_flatten_pattern(::Tuple{}) = []
_flatten_pattern(p::Union{AbstractVector,Tuple}) = Iterators.flatten(_flatten_pattern.(p))
_flatten_pattern(p::Not) = [p.pattern]
_flatten_pattern(p::Pattern) = [p]
