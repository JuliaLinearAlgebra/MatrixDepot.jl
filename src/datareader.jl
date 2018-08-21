# return info comment strings for UF sparse matrix
if VERSION >= v"0.7.0-"
    Sparse = SuiteSparse.CHOLMOD.Sparse
elseif isdefined(:SparseArrays)
    Sparse = Base.SparseArrays.CHOLMOD.Sparse
else
    Sparse = Base.SparseMatrix.CHOLMOD.Sparse
end

function ufinfo(filename::AbstractString)
    io = IOBuffer()
    open(filename,"r") do mmfile
        ll = readline(mmfile)
        while length(ll) > 0 && ll[1] == '%'
            println(io, ll)
            ll = readline(mmfile)
        end
    end
    String(take!(io))
end

"""
    mreader(dat::RemoteMatrixData)
return dictionary with all data (matrices, rhs, other metadata).
The data may be accessed using the key contained in data.metadata.
"""
function mreader(data::RemoteMatrixData)
    mdata = data.datacache.value
    if mdata === nothing
        if length(data.metadata) == 0
            loadmatrix(data)
        end
        mdata = Dict{AbstractString,Any}()
        for name in data.metadata
            path = joinpath(localdir(data), name)
            mdata[name] = endswith(name, ".mtx") ? sparse(Sparse(path)) : read(path, String)
        end
        data.datacache.value = mdata
    end
    mdata
end
mreader(data::MatrixData) = nothing

"""
    mreader(data:RemoteMatrixData, key::AbstractString)
return specific data files (matrix, rhs, solution, or other metadata.
The `key` must be contained in data.metadata or `nothing` is returned.
"""
function mreader(data::RemoteMatrixData, metaname::AbstractString)
    mdata = reader(data)
    mdata !== nothing ? mdata[metaname] : nothing
end

#internal helper to select special metadata (matrix, rhs, or solution)
function readmetaext(data::RemoteMatrixData, exli::AbstractString...)
    base = rsplit(data.name, '/', limit=2)[end]
    mdata = mreader(data)
    for ext in exli
        f = string(base, '_', ext, ".mtx")
        if haskey(mdata, f)
            return mdata[f]
        end
    end
    nothing
end

matrix(data::RemoteMatrixData) = readmetaext(data, "")
rhs(data::RemoteMatrixData) = readmetaext(data, "_b", "_rhs1", "_rhs")
solution(data::RemoteMatrixData) = readmetaext(data, "_x")

matrix(data::GeneratedMatrixData, args...) = data.func(args...)
matrix(data::MatrixData, args...) = nothing
rhs(data::MatrixData, args...) = nothing
solution(data::MatrixData, args...) = nothing

