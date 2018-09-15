
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
    m::Int
    n::Int
    nnz::Int
    dnz::Int
    kind::String
    date::Int # the Julian year number e.g. 2018
    title::AbstractString
    author::AbstractString
    ed::AbstractString
    fields::AbstractString
    notes::AbstractString
end

## Matrix objects

struct RemoteParameters
    dataurl::String
    indexurl::String
    indexgrep::String
    scan::Tuple
    extension::String
end

abstract type RemoteType end

struct TURemoteType <: RemoteType
    params::RemoteParameters
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
end
Base.show(io::IO, db::MatrixDatabase) = print(io, "MatrixDatabase(", length(db.data), ")")

# Patterns

const IntOrVec = Union{Integer,AbstractVector{<:Integer}}

struct Alias{T,D<:Union{IntOrVec,Colon,AbstractVector{<:IntOrVec}}}
    data::D
    Alias{T}(d::D) where {T<:MatrixData,D} = new{T,D}(d)
    Alias{T}(d...) where {T<:MatrixData} = new{T,Vector{<:IntOrVec}}(Vector{IntOrVec}(collect(d)))
end

abstract type AbstractNot end

const Pattern = Union{Function,AbstractString,Regex,Symbol,Alias,
                       AbstractVector,Tuple,AbstractNot}

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
    print(io, " ", data.kind, " ")
    meta = join(metastring.(data.name, metadata(data)), ", ")
    n = length(meta)
    if n > 40
        meta = string(meta[1:17], " ... ", meta[end-17:end])
    end
    print(io, "[", meta, "]")
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
aliasname(::Type{RemoteMatrixData{TURemoteType}}, i::Integer) = string('#', i)
aliasname(::Type{RemoteMatrixData{MMRemoteType}}, i::Integer) = string('#', 'M', i)
aliasname(::Type{GeneratedMatrixData{N}}, i::Integer) where N = string('#', N, i)
function aliasname(T::Type{<:MatrixData}, r::AbstractVector{<:IntOrVec})
    aliasname.(T, [ x for x in Iterators.flatten(r) if x > 0])
end
aliasname(T::Type{<:MatrixData}, r::Colon) = aliasname(T, 0)
aliasname(data::MatrixData) = aliasname(typeof(data), data.id)
aliasname(ali::Alias{T,D}) where {T,D} = aliasname(T, ali.data)

Base.show(io::IO, a::Alias) = print(io, aliasname(a))

import Base: get, empty!
get(db::MatrixDatabase, key::Tuple, default=nothing) = get(db.data, key, default)
function get(db::MatrixDatabase, name::AbstractString, default=nothing)
    get(db.data, keyfor(db, name), default)
end

empty!(db::MatrixDatabase) = (empty!(db.aliases); empty!(db.data))

function keyfor(db::MatrixDatabase, name::AbstractString)
    key = Tuple(split(name, '/'))
    if haskey(db.data, name)
        name
    elseif length(key) == 1
        get(db.aliases, name, nothing)
    else
        name
    end
end

localindex(::TURemoteType) = abspath(DATA_DIR, "uf_matrices.html")
localindex(::MMRemoteType) = abspath(DATA_DIR, "mm_matrices.html")

localbase(::Type{TURemoteType}) = abspath(DATA_DIR, "uf")
localbase(::Type{MMRemoteType}) = abspath(DATA_DIR, "mm")

localdir(data::RemoteMatrixData{T}) where T = abspath(localbase(T), filename(data))
localdir(data::MatrixData) = nothing

indexurl(remote::RemoteType) = remote.params.indexurl
dataurl(remote::RemoteType) = remote.params.dataurl
dataurl(::Type{T}) where T<:RemoteType = dataurl(preferred(T))
extension(::Type{T}) where T<:RemoteType = extension(preferred(T))
extension(remote::RemoteType) = remote.params.extension
dataurl(data::RemoteMatrixData{T}) where T = join((dataurl(T), data.name), '/') * extension(T)
filename(data::RemoteMatrixData) = joinpath(split(data.name, '/')...)
localfile(data::RemoteMatrixData{TURemoteType}) = localdir(data) * ".tar.gz"
localfile(data::RemoteMatrixData{MMRemoteType}) = localdir(data) * ".mtx.gz"
function matrixfile(data::RemoteMatrixData{TURemoteType})
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

#=
row_num(data::RemoteMatrixData) = data.header.m
col_num(data::RemoteMatrixData) = data.header.n
nz_num(data::RemoteMatrixData) = data.header.nnz
dnz_num(data::RemoteMatrixData) = data.header.dnz
kind(data::RemoteMatrixData) = data.header.kind
date(data::RemoteMatrixData) = data.header.date
row_num(data::MatrixData) = 0
col_num(data::MatrixData) = 0
nz_num(data::MatrixData) = 0
dnz_num(data::MatrixData) = 0
kind(data::MatrixData) = ""
date(data::MatrixData) = 0
ident(data::MatrixData) = data.id
=#

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

Base.length(a::Alias) = length(a.data)
function Base.iterate(a::Alias{T}, args...) where T
    d = iterate(a.data, args...); d === nothing && return d
    Alias{T}(d[1]), d[2]
end
