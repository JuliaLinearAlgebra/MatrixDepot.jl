# The implementations are inspired by MatrixMarket.jl
# https://github.com/JuliaSparse/MatrixMarket.jl
# The MatrixMarket.jl package is licensed under the MIT Expat License:
# Copyright (c) 2013: Viral B. Shah.

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#####################################################
# Download data from UF Sparse Matrix Collection
#####################################################

function remotetype(line::AbstractString)
    global uf_remote
    TAMUID = TA_REMOTE.params.indexgrep
    UFID = UF_REMOTE.params.indexgrep
    MMID = MM_REMOTE.params.indexgrep
    if occursin(TAMUID, line)
        uf_remote = TA_REMOTE
        uf_remote
    elseif occursin(UFID, line)
        uf_remote = UF_REMOTE
        uf_remote
    elseif occursin(MMID, line)
        MM_REMOTE
    else
        nothing
    end
end

# extract loading url base and matrix names from index file
function extract_names(matrices::AbstractString, db::MatrixDatabase=MATRIX_DB)
    remote = nothing
    matrixdata = RemoteMatrixData[]
    count = 0
    open(matrices) do f
        id = nothing
        for line in readlines(f)
            if remote === nothing
                remote = remotetype(line)
                continue
            end
            atyp, grepex, spquote, ending, parts, regexid = remote.params.scan
            m = regexid === nothing ? nothing : match(regexid, line)
            if m != nothing
                id = m.captures[1]
            end
            if occursin(grepex, line)
                murl = split(line, '"')[spquote]
                if endswith(murl, ending)
                    list = rsplit(murl[1:end-length(ending)], '/', limit=parts+1)[2:end]
                    count += 1
                    name = join(list, '/')
                    alias = id === nothing ? string(count) : id
                    le = string('%', list[end])
                    db.aliases["#$atyp$alias"] = name
                    while Base.get(db.aliases, le, nothing) !== nothing
                        le = string('%', le)
                    end
                    db.aliases[le] = name
                    count = id != nothing ? parse(Int, id) : count
                    data = RemoteMatrixData(name, le, count, remote)
                    addmetadata!(data)
                    push!(db, data)
                    id = nothing
                end
            end
        end
    end
    matrixdata
end

function downloadindex(remote::RemoteType)
    file = localindex(remote)
    url = indexurl(remote)
    isfile(file) || download(url, file)
    nothing
end

# download html files and store matrix data as a list of tuples
function downloadindices(db::MatrixDatabase; generate_list::Bool = true)
    # UF Sparse matrix collection
    global uf_remote
    empty!(db)
    try
        downloadindex(preferred_uf())
    catch # fallback if first site does not succeed
        uf_remote = alternate_uf()
        downloadindex(uf_remote)
    end

    downloadindex(MM_REMOTE)

    if generate_list
        md = matrix_data_name_list()
        uf_matrixdata = extract_names(localindex(preferred_uf()), db)
        mm_matrixdata = extract_names(localindex(MM_REMOTE), db)
    end
    nothing
end

# given a matrix name, 
function search(matrixname::AbstractString, db::MatrixDatabase=MATRIX_DB)
    uf_matrixdata, mm_matrixdata = downloadindices(db)
    uf_matrices = AbstractString[]
    mm_matrices = AbstractString[]
    [push!(uf_matrices, m[2]) for m in uf_matrixdata]
    [push!(mm_matrices, m[3]) for m in mm_matrixdata]

    datalist = AbstractString[]
    if (matrixname in uf_matrices) || (matrixname in mm_matrices)
        uf_index = findall(isequal(matrixname), uf_matrices)
        mm_index = findall(isequal(matrixname), mm_matrices)
        for i in uf_index
            collectionname, matrixname = uf_matrixdata[i]
            push!(datalist, string(collectionname, '/', matrixname))
        end
        for i in mm_index
            collectionname, setname, matrixname = mm_matrixdata[i]
            push!(datalist, string(collectionname, '/', setname, '/', matrixname))
        end
        return datalist
    else
        error("can not find $matrixname in the database.")
    end
end

function gunzip(fname)
    endswith(fname, ".gz") || error("gunzip: $fname: unknown suffix")

    destname = rsplit(fname, ".gz", limit=2)[1]

    open(destname, "w") do f
        GZip.open(fname) do g
            write(f, read(g))
        end
    end
    destname
end

function matchnames(p, a)
    n = length(a)
    length(p) == n || (return false)
    for i = 1:n
        p[i] in ("*", a[i]) || (return false)
    end
    true
end

# loadmatrix
# --------------
# loadmatrix(NAME) download a matrix from UF or MM sparse matrix collection
# where NAME is a string of collection name + '/' + matrix name.
#
# Example
# -------
# MatrixDepot.loadmatrix("HB/1138_bus") # uf sparse matrix
# MatrixDepot.loadmatrix("Harwell-Boeing/psadmit/1138_bus") # matrix market
#
function loadmatrix(name::AbstractString, db::MatrixDatabase=MATRIX_DB)   
    data = get(db, name)
    file = matrixfile(data)
    if isfile(file)
        return
    end
    dirfn = localfile(data)
    dir = dirname(localdir(data))
    url = dataurl(data)

    isdir(dir) || mkpath(dir)

    try
        println("downloading: ", url)
        download(url, dirfn)
        tarfile = gunzip(dirfn)
        if endswith(tarfile, ".tar")
            run(`tar -vxf $tarfile -C $dir`)
            rm(tarfile)
        end
    finally
        rm(dirfn)
    end
    nothing
end

function loadmatrix(data::MatrixData, db::MatrixDatabase=MATRIX_DB)   
    stringvec = search(name)
    if length(stringvec) == 1
        return matrixdepot(stringvec[1], :get)
    else
        println("Try loadmatrix(`name`), where `name` is one of the elements in the following Array:")
        return stringvec
    end
end

function addmetadata!(data::RemoteMatrixData)
    file = matrixfile(data)
    dir = dirname(file)
    empty!(data.metadata)
    isdir(dir) || return
    base = basename(file)
    name, ext = rsplit(base, '.', limit=2)
    filtop(x) = x == base || startswith(x, string(name, '_'))
    append!(data.metadata, filter(filtop, readdir(dir)))
    nothing
end
