
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

struct RemoteMatrixData <:MatrixData
    name::AbstractString
    alias::String
    id::Int
    metadata::Vector{AbstractString}
    remote::RemoteType
    infocache::WeakRef
    datacache::WeakRef
    RemoteMatrixData(name, alias, id, rem) =
    new(name, alias, id, AbstractString[], rem, WeakRef(), WeakRef())
end

struct MatrixDatabase
    data::Dict{AbstractString,MatrixData}
    aliases::Dict{AbstractString,AbstractString}
    MatrixDatabase() = new(Dict{AbstractString,MatrixData}(),
                           Dict{AbstractString,AbstractString}())
end

function Base.push!(db::MatrixDatabase, data::MatrixData)
    key = data.name
    alias = data.alias
    db.data[data.name] = data
    db.aliases[data.alias] = key
    db.aliases[aliasid(data)] = key
end

aliasid(data::RemoteMatrixData) = aliasid(data.remote, data.id)
aliasid(::TURemoteType, id) = string('#', id)
aliasid(::MMRemoteType, id) = string('#', 'M', id)

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

localbase(::TURemoteType) = abspath(DATA_DIR, "uf")
localbase(::MMRemoteType) = abspath(DATA_DIR, "mm")

localdir(data::RemoteMatrixData) = abspath(localbase(data.remote), filename(data))
localdir(data::MatrixData) = nothing

indexurl(remote::RemoteType) = remote.params.indexurl
dataurl(remote::RemoteType) = remote.params.dataurl
extension(remote::RemoteType) = remote.params.extension
dataurl(data::RemoteMatrixData) = join((dataurl(data.remote), data.name), '/') * extension(data.remote)
filename(data::RemoteMatrixData) = joinpath(split(data.name, '/')...)
localfile(data::RemoteMatrixData) = localfile(data, data.remote)
localfile(data::RemoteMatrixData, ::TURemoteType) = localdir(data) * ".tar.gz"
localfile(data::RemoteMatrixData, ::MMRemoteType) = localdir(data) * ".mtx.gz"
matrixfile(data::RemoteMatrixData) = matrixfile(data, data.remote)
matrixfile(data::RemoteMatrixData, ::TURemoteType) = joinpath(localdir(data), rsplit(data.name, '/', limit=2)[end] * ".mtx")
matrixfile(data::RemoteMatrixData, ::MMRemoteType) = localdir(data) * ".mtx"
