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

# check line for remote identification and store in global variable
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

function parse_headerinfo(akku::Dict{AbstractString,AbstractString}, count::Integer)
    toint(id::String) = parse(Int, replace(get(akku, id, "0"), ','=>""))
    id = toint("id")
    id == 0 && (id = count)
    m = toint("num_rows")
    n = toint("num_cols")
    nnz = toint("nonzeros")
    kind = get(akku, "kind", "")
    datestr = get(akku, "date", "0")
    date = match(r"^\d+$", datestr) != nothing ? parse(Int, datestr) : 0
    id, IndexInfo(m, n, nnz, 0, kind, date)
end

# extract loading url base and matrix names from index file
function extract_names(db::MatrixDatabase, matrices::AbstractString)
    remote = nothing
    matrixdata = RemoteMatrixData[]
    count = 0
    open(matrices) do f
        akku = Dict{AbstractString,AbstractString}()
        for line in readlines(f)
            if remote === nothing
                remote = remotetype(line)
                continue
            end
            _, grepex, spquote, ending, parts, regexinf = remote.params.scan
            m = regexinf === nothing ? nothing : match(regexinf, line)
            if m != nothing
                pname = length(m.captures) == 2 ? m.captures[1] : "id"
                akku[pname] = m.captures[end]
            end
            if occursin(grepex, line)
                murl = split(line, '"')[spquote]
                if endswith(murl, ending)
                    list = rsplit(murl[1:end-length(ending)], '/', limit=parts+1)[2:end]
                    count += 1
                    name = join(list, '/')
                    id, info = parse_headerinfo(akku, count)
                    data = RemoteMatrixData{typeof(remote)}(name, id, info)
                    addmetadata!(data)
                    push!(db, data)
                    count = id
                end
            end
        end
    end
    matrixdata
end

# read indexfile
function downloadindex(remote::RemoteType)
    file = localindex(remote)
    url = indexurl(remote)
    if !isfile(file)
        println("dowloading index file $url")
        download(url, file)
    end
    nothing
end

# collect the keys from local database (MATRIXDICT od USERMATRIXDICT)
# provide a numerical id counting from 1 for either database.
function insertlocal(db::MatrixDatabase, T::Type{<:GeneratedMatrixData},
                     ldb::Dict{<:AbstractString,Function})

    cnt = 0
    ks = sort!(collect(keys(ldb)))
    for k in ks
        cnt += 1
        db.data[k] = T(k, cnt, ldb[k])
    end
end

"""
    MatrixDepot.downloadindices(db)
download html files and store matrix data in `db`.
additionally generate aliases for local and loaded matrices.
"""
function downloadindices(db::MatrixDatabase)
    # UF Sparse matrix collection
    global uf_remote
    empty!(db)
    insertlocal(db, GeneratedMatrixData{:B}, MATRIXDICT)
    insertlocal(db, GeneratedMatrixData{:U}, USERMATRIXDICT)

    try
        downloadindex(preferred(TURemoteType))
        downloadindex(preferred(MMRemoteType))

        extract_names(db, localindex(preferred(TURemoteType)))
        extract_names(db, localindex(preferred(MMRemoteType)))
    finally
        for data in values(db.data)
            db.aliases[aliasname(data)] = data.name
        end
    end
    nothing
end

"""
    MatrixDepot.update()
update database index from the websites
"""
function update(db::MatrixDatabase=MATRIX_DB)
    uf_matrices = localindex(preferred(TURemoteType))
    isfile(uf_matrices) && rm(uf_matrices)
    mm_matrices = localindex(preferred(MMRemoteType))
    isfile(mm_matrices) && rm(mm_matrices)
    downloadindices(db)
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
function loadmatrix(db::MatrixDatabase, data::RemoteMatrixData)
    file = matrixfile(data)
    if isfile(file)
        return 0
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
    1
end
loadmatrix(db::MatrixDatabase, data::GeneratedMatrixData) = 0

function data_warn(data::RemoteMatrixData, dn, i1, i2)
    @warn "$(data.name): header $dn = $i1 file $dn = $i2 $(data.properties[])"
    i1
end

function extremnnz(data::RemoteMatrixData, m, n, k)
    if issymmetric(data) || ishermitian(data) || isskew(data)
        if m != n
            @warn "$(data.name) needs to be quadratic but is ($mx$n)"
        end
        isskew(data) ? (2k, 2k) : (2k - max(m, n), 2k)
    else
        (k, k)
    end
end
function extremnnz(data::RemoteMatrixData)
    extremnnz(data, row_num(data), col_num(data), dnz_num(data))
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

    if (finfo = mmreadheader(file)) !== nothing
        data.properties[] = MMProperties(finfo[4:end]...)
        m, n, k = finfo[1:3]
        n1, n2 = extremnnz(data, m, n, k)
        hdr = data.header
        hdr.m = hdr.m in (0, m) ? m : data_warn(data, "m", hdr.m, m)
        hdr.n = hdr.n in (0, n) ? n : data_warn(data, "n", hdr.n, n)
        hdr.nnz = hdr.nnz == 0  ? n2 : hdr.nnz <= n2 ? hdr.nnz :
                                        data_warn(data, "nnz", hdr.nnz, k)
        hdr.dnz = k
    end
    nothing
end

