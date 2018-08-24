
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
    metadata::Vector{AbstractString}
    infocache::WeakKeyDict{String}
    datacache::WeakKeyDict{String}
    RemoteMatrixData{T}(name, id) where T =
        new(name, id, AbstractString[], WeakKeyDict{String,Any}(), WeakKeyDict{String,Any}())
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
