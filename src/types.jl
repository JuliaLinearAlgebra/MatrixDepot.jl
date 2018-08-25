
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

mutable struct MatrixReference
    name::String
    data::Any
end

abstract type MatrixData end

struct RemoteMatrixData{T<:RemoteType} <:MatrixData
    name::AbstractString
    id::Int
    properties::Ref{Union{<:MMProperties,Nothing}}
    metadata::Vector{AbstractString}
    infocache::Dict{String}
    datacache::Dict{String}
    function RemoteMatrixData{T}(name, id) where T
        new(name, id, Ref{Union{<:MMProperties,Nothing}}(nothing), AbstractString[],
            Dict{String,Any}(), Dict{String,Any}())
    end
end

function Base.show(io::IO, data::RemoteMatrixData)
    print(io, "RemoteMatrixData(")
          print(io, "$(data.name)($(data.id)) ")
    show(io, data.properties[])
    print(io, " [")
    print(io, join(data.metadata, ", "))
    print(io,"])")
end

abstract type GeneratedMatrixData <:MatrixData end

struct GeneratedBuiltinMatrixData <:GeneratedMatrixData
    name::AbstractString
    id::Int
    func::Function
end 
struct GeneratedUserMatrixData <:GeneratedMatrixData
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

function Base.push!(db::MatrixDatabase, data::MatrixData)
    key = data.name
    db.data[data.name] = data
    db.aliases[aliasid(data)] = key
end

aliasid(i::Integer) = string('#', i)
aliasid(data::RemoteMatrixData{TURemoteType}) = aliasid(data.id)
aliasid(data::RemoteMatrixData{MMRemoteType}) = string('#', 'M', data.id)

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
matrixfile(data::RemoteMatrixData{TURemoteType}) =
    joinpath(localdir(data), rsplit(data.name, '/', limit=2)[end] * ".mtx")
matrixfile(data::RemoteMatrixData{MMRemoteType}) = localdir(data) * ".mtx"

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
    print(io, "{")
    show(io, mp.object); print(io, ", ")
    show(io, mp.format); print(io, ", ")
    show(io, mp.field); print(io, ", ")
    show(io, mp.symmetry); print(io, "}")
end


