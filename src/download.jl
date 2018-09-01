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
                    db.aliases["#$atyp$alias"] = name
                    count = id != nothing ? parse(Int, id) : count
                    data = RemoteMatrixData{typeof(remote)}(name, count)
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
    if !isfile(file)
        println("dowloading index file $url")
        download(url, file)
    end
    nothing
end

function insertlocal(T::Type{<:GeneratedMatrixData},
                     ldb::Dict{<:AbstractString,Function},
                     db::MatrixDatabase)

    cnt = 0
    ks = sort!(collect(keys(ldb)))
    for k in ks
        cnt += 1
        db.data[k] = T(k, cnt, ldb[k])
    end
end

# download html files and store matrix data in MATRIX_DB
function downloadindices(db::MatrixDatabase)
    # UF Sparse matrix collection
    global uf_remote
    empty!(db)
    insertlocal(GeneratedBuiltinMatrixData, MATRIXDICT, db)
    insertlocal(GeneratedUserMatrixData, USERMATRIXDICT, db)
    try
        downloadindex(preferred(TURemoteType))
    catch # fallback if first site does not succeed
        uf_remote = alternate(TURemoteType)
        downloadindex(uf_remote)
    end

    downloadindex(preferred(MMRemoteType))

    extract_names(localindex(preferred(TURemoteType)), db)
    extract_names(localindex(preferred(MMRemoteType)), db)
    nothing
end

# update database from the websites
function update()
    uf_matrices = string(DATA_DIR, "/uf_matrices.html")
    mm_matrices = string(DATA_DIR, "/mm_matrices.html")
    if isfile(uf_matrices)
        rm(uf_matrices)
    end
    if isfile(mm_matrices)
        rm(mm_matrices)
    end
    downloadindices(MATRIX_DB)
end

function gunzip(fname)
    endswith(fname, ".gz") || error("gunzip: $fname: unknown suffix")

    destname = rsplit(fname, ".gz", limit=2)[1]
    BUFFSIZE = 1000000
    open(destname, "w") do f
        GZip.open(fname) do g
            buffer = read(g, BUFFSIZE)
            while length(buffer) > 0
                write(f, buffer)
                buffer = read(g, BUFFSIZE)
            end
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
function loadmatrix(data::RemoteMatrixData, db::MatrixDatabase)   
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
            rm(tarfile; force=true)
        end
    finally
        rm(dirfn, force=true)
    end
    addmetadata!(data)
    nothing
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

    if (finfo = fileinfo(file)) !== nothing
        data.properties[] = finfo
    end
    nothing
end

