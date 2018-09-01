# return info comment strings for UF sparse matrix

function ufinfo(filename::AbstractString)
    io = IOBuffer()
    open(filename,"r") do mmfile
        line = readline(mmfile)
        while startswith(line, '%')
            println(io, line)
            line = readline(mmfile)
        end
        while isempty(strip(line))
            line = readline(mmfile)
        end
        println(io, line)
    end
    String(take!(io))
end

function trymmread(path::AbstractString)
    mmread(path)
end

"""
    metareader(dat::RemoteMatrixData)
return dictionary with all data (matrices, rhs, other metadata).
The data may be accessed using the keys contained in data.metadata.
The keys are identical to the names of the files keeping the data.
"""
function metareader(data::RemoteMatrixData)
    result = Dict{String,Any}()
    for name in data.metadata
        result[name] = metareader(data, name)
    end
    result
end

metareader(data::MatrixData) = nothing

"""
    metareader(data:RemoteMatrixData, key::AbstractString)
return specific data files (matrix, rhs, solution, or other metadata.
The `key` must be contained in data.metadata or `nothing` is returned.
"""
function metareader(data::RemoteMatrixData, name::AbstractString)
    dc = data.datacache
    if haskey(dc, name)
        return dc[name]
    end
    if name in data.metadata
        path = joinpath(dirname(matrixfile(data)), name)
        m = endswith(name, ".mtx") ? trymmread(path) : read(path, String)
        if data.status[]
            dc[name] = m
        end
        m
    else
        nothing
    end
end

#internal helper to select special metadata (matrix, rhs, or solution)
function readmetaext(data::RemoteMatrixData, exli::AbstractString...)
    base = rsplit(data.name, '/', limit=2)[end]
    for ext in exli
        f = string(base, ext, ".mtx")
        if f in data.metadata
            return metareader(data, f)
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

