"""
    metareader(data:RemoteMatrixData, key::AbstractString)
return specific data files (matrix, rhs, solution, or other metadata.
The `key` must be contained in data.metadata or an `DataError` is thrown.
"""
function metareader(data::RemoteMatrixData, name::AbstractString)
    if name in data.metadata
        path = joinpath(dirname(matrixfile(data)), name)
        endswith(name, ".mtx") ? mmread(path) : read(path, String)
    else
        daterr("$(data.name) has no metadata $name")
    end
end

function metareader(mdesc::MatrixDescriptor, name::AbstractString)
    get!(mdesc.cache, name) do
        metareader(mdesc.data, name)
    end
end

#internal helper to select special metadata (matrix, rhs, or solution)
function metaname(data::RemoteMatrixData, exli::AbstractString...)
    base = rsplit(data.name, '/', limit=2)[end]
    for ext in exli
        f = string(base, ext, ".mtx")
        if f in data.metadata
            return f
        end
    end
    daterr("unknown metadata extensions: `$([exli...])` - available $(String.(data.metadata))")
end

function readmetaext(data::RemoteMatrixData, exli::AbstractString...)
    metareader(data, metaname(data, exli...))
end

matrix(data::RemoteMatrixData) = readmetaext(data, "")
rhs(data::RemoteMatrixData) = readmetaext(data, "_b", "_rhs1", "_rhs")
solution(data::RemoteMatrixData) = readmetaext(data, "_x")

matrix(data::GeneratedMatrixData, args...) = data.func(args...)
matrix(data::MatrixData, args...) = throw(MethodError(matrix, (data, args...)))
rhs(data::MatrixData, args...) = throw(MethodError(matrix, (data, args...)))
solution(data::MatrixData, args...) = throw(MethodError(matrix, (data, args...)))

function fillcache!(mdesc::MatrixDescriptor{<:GeneratedMatrixData})
    x = mdesc.cache[]
    if x === nothing
        x = mdesc.cache[] = mdesc.data.func(mdesc.args...)
    end
    x
end

matrix(mdesc::MatrixDescriptor{<:GeneratedMatrixData}) = get_A(fillcache!(mdesc))
rhs(mdesc::MatrixDescriptor{<:GeneratedMatrixData}) = get_b(fillcache!(mdesc))
solution(mdesc::MatrixDescriptor{<:GeneratedMatrixData}) = get_x(fillcache!(mdesc))

get_A(a::AbstractMatrix) = a
get_A(a::Union{RegProb,RegProbNoSolution}) = a.A
get_b(a::Union{RegProb,RegProbNoSolution}) = a.b
get_x(a::RegProb) = a.x
get_A(a) = argerr("no matrix available")
get_b(a) = argerr("no rhs available")
get_x(a) = argerr("no solution available")

matrix(mdesc::MatrixDescriptor{<:RemoteMatrixData}) = readmetaext(mdesc, "")
rhs(mdesc::MatrixDescriptor{<:RemoteMatrixData}) = readmetaext(mdesc, "_b", "_rhs", "_rhs1")
solution(mdesc::MatrixDescriptor{<:RemoteMatrixData}) = readmetaext(mdesc, "_x")

function readmetaext(mdesc::MatrixDescriptor, exli::AbstractString...)
    f = metaname(mdesc.data, exli...)
    get!(mdesc.cache, f) do
        readmetaext(mdesc.data, exli...)
    end
end

function metareader(mdesc::MatrixDescriptor{<:RemoteMatrixData}, name::AbstractString)
    get!(mdesc.cache, name) do
        metareader(mdesc.data, name)
    end
end

import Base.getproperty
function getproperty(mdesc::MatrixDescriptor{T}, s::Symbol) where T
    s in (:data, :args, :cache) && return getfield(mdesc, s)
    if s == :A
        matrix(mdesc)
    elseif s == :b
        rhs(mdesc)
    elseif s == :x
        solution(mdesc)
    elseif T <: RemoteMatrixData
        data = mdesc.data
        f = metaname(data, string('_', s))
        metareader(mdesc, f)
    else
        data = mdesc.data
        try
            getproperty(data, s)
        catch
            argerr("this instance of $(typeof(data)) has no field $s")
        end
    end
end

function meta_file_to_symbol(name::AbstractString, metanam::AbstractString)
    name = split(name, '/')[end]
    meta, metaext = rsplit(metanam, '.', limit=2)
    afteru(meta, name) = meta[nextind(meta, sizeof(name)):end]
    startswith(meta, name) || daterr("$metanam is not a metaname for $name")
    meta = afteru(meta, name)
    ismtx = metaext == ".mtx"
    s2(meta) = meta[2:end]
    s2(meta, ext) = string(s2(meta), '_', ext, '_')
    if meta == ""
        ismtx ? :A : Symbol(metaext, '_')
    elseif match(r"^_[_[:alpha]]", meta) !== nothing
        ismtx ? Symbol(s2(meta)) : Symbol(s2(meta, metaext))
    elseif match(r"^[_[^:alpha]]", meta) !== nothing
        ismtx ? Symbol("_#", s2(meta)) : Symbol("_#", s2(meta, metaext))
    end
end

function meta_symbol_to_file(name::AbstractString, sym::Symbol)
    ssym = string(sym)
    if endswith(ssym, '_')
        adds, ext = rsplit(ssym, '_', limit=3)
        if startswith(adds, "_#")
            adds = add[3:end]
        end
        name = split(name, '/')[end]
        string(name, '_', adds, '.', ext)
    else
        adds = ssym
        if startswith(adds, "_")
            adds = adds[3:end]
        end
        name = split(name, '/')[end]
        string(name, '_', adds, ".mtx")
    end
end


