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
    id, MetaInfo(m, n, nnz, 0, kind, date, "", "", "", "", "")
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
Serialize db object in file `db.data`.
If a file `db.data` is present in the data directory, deserialize
the data instead of downloading and collection data.
"""
function downloadindices(db::MatrixDatabase; ignoredb=false)
    # UF Sparse matrix collection
    cachedb = abspath(DATA_DIR, "db.data")
    if !ignoredb && isfile(cachedb)
        open(cachedb, "r") do io
            dbx = deserialize(io)
            merge!(db.data, dbx.data)
            merge!(db.aliases, dbx.aliases)
        end
    else
        _downloadindices(db)
        open(cachedb, "w") do io
            serialize(io, db)
        end
    end
    nothing
end

function _downloadindices(db::MatrixDatabase)
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

"""
    loadmatrix(data::RemoteMatrixData)

Download the files backing the data from a remote repository. That is currently
the TAMU site for the UFl collection and the NIST site for the MatrixMarket
caoolection. The files are uncompressed and un-tar-ed if necessary.
The data files containing the matrix data have to be in MatrixMarket format in
both cases. Note, that some of the files of the MM collection are not available 
in MatrixMarket format. An error message results, if trie to load them.
"""
function loadmatrix(data::RemoteMatrixData)
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
loadmatrix(data::GeneratedMatrixData) = 0

"""
    loadinfo(data::RemoteDate)
Download the first part of the data file. Stop reading, as soon as the initial
comment an the size values of the main matrix have been finished. Stor this in
a file with extension `.info` in the same directory, where the `.mtx` file is.
If the complete file is already availble, the download is not performed, because
the head of the `.mtx` file contains the same lines.
"""
function loadinfo(data::RemoteMatrixData)
    filemtx = matrixfile(data)
    file = matrixinfofile(data)
    if isfile(filemtx) || isfile(file)
        return 0
    end
    url = dataurl(data)
    urls = rsplit(url, '.', limit=3)
    cmd = []
    push!(cmd, `sh -c 'curl "'$url'" -o - 2>/dev/null'`)
    if urls[end] == "gz"
        push!(cmd, `gunzip`)
        resize!(urls, length(urls)-1)
    end
    if url[end] == "tar"
        push!(cmd, `tar -xOf -`)
    end

    out = IOBuffer()
    s = try
        println("downloading head of $url")
        open(pipeline(cmd...), "r") do io
            while ( s = readline(io) ) != ""
                println(out, s)
                s[1] != '%' && break
            end;
            close(io)
        end
        String(take!(out))
    catch ex
        ex isa InterruptException && rethrow()
        String(take!(out))
    finally
        close(out)
    end
    if !isempty(s)
        mkpath(dirname(file))
        write(file, s)
        1
    else
        0
    end
end
loadinfo(data::MatrixData) = 0

function data_warn(data::RemoteMatrixData, dn, i1, i2)
    @warn "$(data.name): header $dn = '$i1' file $dn = '$i2' $(data.properties[])"
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
    extremnnz(data, data.m, data.n, data.dnz)
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

    if !isfile(file)
        file = matrixinfofile(data)
    end
    if (finfo = mmreadheader(file)) !== nothing
        data.properties[] = MMProperties(MATRIX, finfo[:format], finfo[:field], finfo[:symmetry])
        m, n, k = finfo[:m], finfo[:n], finfo[:nz]
        n1, n2 = extremnnz(data, m, n, k)
        hdr = data.header
        hdr.m = hdr.m in (0, m) ? m : data_warn(data, "m", hdr.m, m)
        hdr.n = hdr.n in (0, n) ? n : data_warn(data, "n", hdr.n, n)
        hdr.nnz = hdr.nnz == 0  ? n2 : hdr.nnz <= n2 ? hdr.nnz :
                                        data_warn(data, "nnz", hdr.nnz, k)
        hdr.dnz = k
        date = get(finfo, :date, 0)
        hdr.date = hdr.date in (0, date) ? date : data_warn(data, "date", hdr.date, date)
        kind = get(finfo, :kind, "")
        if kind != ""
            if hdr.kind == ""
                hdr.kind = kind
            else
                if lowercase(replace(kind, '-' => ' ')) != lowercase(hdr.kind)
                    data_warn(data, "kind", hdr.kind, kind)
                end
            end
        end
        for s in (:title, :author, :ed, :fields, :notes)
            val = get(finfo, s, "")
            val != "" && setfield!(hdr, s, val)
        end
    end
    nothing
end

