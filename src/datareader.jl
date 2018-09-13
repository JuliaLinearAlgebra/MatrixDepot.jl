"""
    metareader(md::MatrixDescriptor, key::AbstractString)

Return specific data files (matrix, rhs, solution, or other metadata).
The `key` must be contained in data.metadata or an `DataError` is thrown.
Hint: `md.A, md.b, md.x` also deliver the matrix, rhs, and solution.
So this is needed specifically for other metadata.
"""
function metareader(mdesc::MatrixDescriptor{<:RemoteMatrixData}, name::AbstractString)
    name = metastring_reverse(mdesc.data, name)
    get!(mdesc.cache, name) do
        _metareader(mdesc.data, name)
    end
end

function metareader(mdesc::MatrixDescriptor{<:GeneratedMatrixData}, name::AbstractString)
    fillcache!(mdesc)
    dao = mdesc.cache[]
    if name == "A" && dao isa AbstractArray
        dao
    else
        try
            getproperty(dao, Symbol(name))
        catch
            daterr("this instance of '$(typeof(dao))' has no property '$name'")
        end
    end
end

function _metareader(data::RemoteMatrixData, name::AbstractString)
    if name in data.metadata
        path = joinpath(dirname(matrixfile(data)), name)
        endswith(name, ".mtx") ? mmread(path) : read(path, String)
    else
        daterr("$(data.name) has no metadata $name")
    end
end

function fillcache!(mdesc::MatrixDescriptor{<:GeneratedMatrixData})
    dao = mdesc.cache[]
    if dao === nothing
        dao = mdesc.cache[] = mdesc.data.func(mdesc.args...)
        dao !== nothing || daterr("function $(mdesc.data.func) returned `nothing`")
    end
    dao
end

# This a the preferred API to access metadata.
import Base.getproperty

function getproperty(mdesc::MatrixDescriptor, s::Symbol)
    s in (:data, :args, :cache) && return getfield(mdesc, s)
    metareader(mdesc, string(s))
end

#internal helper to select special metadata (matrix, rhs, or solution)
function metaname(data::RemoteMatrixData, exli::AbstractString...)
    base = rsplit(data.name, '/', limit=2)[end]
    meda = metadata(data)
    for ext in exli
        f = metastring_reverse(data, ext)
        if f in meda
            return f
        end
    end
    daterr("unknown metadata extensions: $exli - available $(String.(data.metadata))")
end

"""
    metastring(name, metaname)

In the standard cases convert full meta name to abbreviated form.

+ given `name = "abc"`, then
+ `"abc.mtx" => "A"`
+ `"abc_def.mtx" => "def"`
+ `"abc_def.txt" => "def.txt"`
+ `"abc-def.mtx" => "abc-def.mtx"`
"""
function metastring(name::AbstractString, metaname::AbstractString)
    MTX = ".mtx"
    base = split(name, '/')[end]
    meta = metaname
    if meta == string(base, MTX)
        meta = "A"
    elseif startswith(meta, string(base, '_'))
        meta = meta[sizeof(base)+2:end]
    end
    es = rsplit(meta, '.', limit=2)
    if endswith(meta, MTX) && meta != metaname
        meta = meta[1:end-sizeof(MTX)]
    end
    meta
end

"""
    metastring_reverse(data::MatrixData, abbreviation::String)

Given a `MatrixData` object and an abbreviation, return a matchingfull metadata name
for the abbreviation. If no full name is found, return the abbreviation unchanged.
"""
function metastring_reverse(data::RemoteMatrixData, metaabbr::AbstractString)
    MTX = ".mtx"
    base = split(data.name, '/')[end]
    metaabbr == "A" && return string(base, MTX)
    mdata = metadata(data)
    meta = string(base, '_', metaabbr)
    meta in mdata && return meta
    meta = string(meta, MTX)
    meta in mdata && return meta
    metaabbr
end
metastring_reverse(::MatrixData, metaabbr::AbstractString) = metaabbr

