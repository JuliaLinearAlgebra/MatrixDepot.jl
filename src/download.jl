#####################################################
# Download data from UF Sparse Matrix Collection
#####################################################

using Base.Filesystem
using ChannelBuffers

# collect the keys from local database (MATRIXDICT or USERMATRIXDICT)
# provide a numerical id counting from 1 for either database.
function insertlocal(db::MatrixDatabase, T::Type{<:GeneratedMatrixData},
                     ldb::Dict{<:AbstractString,Function})

    cnt = 0
    ks = sort!(collect(keys(ldb)))
    for k in ks
        cnt += 1
        push!(db, T(k, cnt, ldb[k]))
    end
    cnt
end

dbpath(db::MatrixDatabase) = abspath(data_dir(), "db.data")
function readdb(db::MatrixDatabase)
    @info("reading database")
    cachedb = dbpath(db)
    open(cachedb, "r") do io
        dbx = deserialize(io)
        merge!(db.data, dbx.data)
        merge!(db.aliases, dbx.aliases)
    end
end

function writedb(db::MatrixDatabase)
    cachedb = dbpath(db)
    #Avoid serializing user matrices, since we don't know which ones will be loaded on deserialization.
    dbx_data = filter(((key, val),) -> !(val isa GeneratedMatrixData{:U}), db.data)
    dbx_aliases = filter(((_, key),) -> !(db.data[key] isa GeneratedMatrixData{:U}), db.aliases)
    dbx = MatrixDatabase(dbx_data, dbx_aliases)
    open(cachedb, "w") do io
        serialize(io, dbx)
    end
end

"""
    MatrixDepot.downloadindices(db)
download index files and store matrix data in `db`.
additionally generate aliases for local and loaded matrices.
Serialize db object in file `db.data`.
If a file `db.data` is present in the data directory, deserialize
the data instead of downloading and collecting data.
"""
function downloadindices(db::MatrixDatabase; ignoredb=false)
    added = 0
    if !ignoredb && isfile(dbpath(db))
        try
            readdb(db)
        catch ex
            println(ex)
            @warn("recreating database file")
            _downloadindices(db)
            added = 1
        end
    else
        @info("creating database file")
        _downloadindices(db)
        added = 1
    end
    @info("adding metadata...")
    added += addmetadata!(db)
    @info("adding svd data...")
    added += addsvd!(db)
    added += insertlocal(db, GeneratedMatrixData{:B}, MATRIXDICT)
    added += insertlocal(db, GeneratedMatrixData{:U}, USERMATRIXDICT)
    if added > 0
        @info("writing database")
        writedb(db)
    end
    nothing
end

function _downloadindices(db::MatrixDatabase)
    @info("reading index files")
    empty!(db)

    try
        readindex(preferred(SSRemoteType), db)
        readindex(preferred(MMRemoteType), db)
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
    cachedb = abspath(data_dir(), "db.data")
    isfile(cachedb) && rm(cachedb)
    uf_matrices = localindex(preferred(SSRemoteType))
    isfile(uf_matrices) && rm(uf_matrices)
    mm_matrices = localindex(preferred(MMRemoteType))
    isfile(mm_matrices) && rm(mm_matrices)
    downloadindices(db, ignoredb=true)
end

"""
    loadmatrix(data::RemoteMatrixData)

Download the files backing the data from a remote repository. That is currently
the TAMU site for the UFl collection and the NIST site for the MatrixMarket
collection. The files are uncompressed and un-tar-ed if necessary.
The data files containing the matrix data have to be in MatrixMarket format in
both cases. Note, that some of the files of the MM collection are not available
in MatrixMarket format. An error message results, if tried to load them.
"""
function loadmatrix(data::RemoteMatrixData)
    file = matrixfile(data)
    if isfile(file)
        addmetadata!(data)
        return 0
    end
    dir = dirname(localdir(data))
    dirt = string(dir, ".tmp")
    rm(dirt, force=true, recursive=true)
    mkpath(dirt)
    url = redirect(dataurl(data))
    isdir(dir) || mkpath(dir)
    wdir = pwd()
    try
        @info("downloading: $url")
        pipe = downloadpipeline(url, dirt)
        wait(run(pipe))
        cd(dirt)
        for file in readdir(dirt)
            mv(file, joinpath(dir, file), force=true)
        end
    catch ex
        @warn("download of $url failed: $ex")
    finally
        rm(dirt, force=true, recursive=true)
        cd(wdir)
    end
    addmetadata!(data)
    1
end
loadmatrix(data::GeneratedMatrixData) = 0

"""
    loadinfo(data::RemoteDate)
Download the first part of the data file. Stop reading, as soon as the initial
comment and the size values of the main matrix have been finished. Store this in
a file with extension `.info` in the same directory, where the `.mtx` file is.
If the complete file is already available, the download is not performed, because
the head of the `.mtx` file contains the same lines.
"""
function loadinfo(data::RemoteMatrixData)
    filemtx = matrixfile(data)
    file = matrixinfofile(data)
    if isfile(filemtx) || isfile(file)
        return 0
    end
    url = redirect(dataurl(data))
    io = ChannelPipe()
    pipe = downloadpipeline(url, io)
    tl = run(pipe)
    out = IOBuffer()
    s = try
        @info("downloading head of $url")
        skip = 0
        while ( s = readskip(io) ) != ""
            skip = s[1] == '%'  || isempty(strip(s)) ? 0 : skip + 1
            skip <= 1 && println(out, s)
            if skip == 1 && length(split(s)) == 3
                break
            end
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
        addmetadata!(data)
        1
    else
        0
    end
end
loadinfo(data::MatrixData) = 0

"""
    downloadpipeline(url)

Set up a command pipeline (external processes to download and expand data)
"""
function downloadpipeline(url::AbstractString, dir=nothing)
    urls = rsplit(url, '.', limit=3)
    cmd = Any[ChannelBuffers.curl(url)]
    if urls[end] == "gz"
        push!(cmd, gunzip())
        resize!(urls, length(urls)-1)
        url = url[1:end-3]
    end
    if dir isa AbstractString
        if urls[end] == "tar"
            push!(cmd, tarx(dir))
        else
            file = rsplit(url, '/', limit=2)
            push!(cmd, joinpath(dir, file[end]))
        end
    elseif dir isa ChannelPipe
        push!(cmd, dir)
    end
    pipeline(cmd...)
end

function downloadcommand(url::AbstractString, filename::AbstractString="-")
    curl(url) → (filename == "-" ? stdout : filename)
end

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

function addmetadata!(db::MatrixDatabase)
    added = 0
    for data in values(db.data)
        if isunloaded(data)
            addmetadata!(data)
            added += isloaded(data)
        end
    end
    added
end

addmetadata!(::MatrixData) = 0
function addmetadata!(data::RemoteMatrixData)
    file = matrixfile(data)
    dir = dirname(file)
    empty!(data.metadata)
    isdir(dir) || return 0
    base = basename(file)
    name, ext = rsplit(base, '.', limit=2)
    filtop(x) = x == base || (startswith(x, string(name, '_')) && !endswith(x, "SVD.mat"))
    append!(data.metadata, filter(filtop, readdir(dir)))

    if !isfile(file)
        file = matrixinfofile(data)
    end
    if (finfo = mmreadheader(file)) !== nothing
        data.properties[] = MMProperties(MATRIX, finfo[:format], finfo[:field], finfo[:symmetry])
        m, n, k = finfo[:m], finfo[:n], get(finfo, :nz, 0)
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
            if hdr.kind == "" || hdr.kind == "other problem"
                hdr.kind = kind
            else
                if standardkind(kind) != standardkind(hdr.kind)
                    data_warn(data, "kind", hdr.kind, kind)
                end
            end
        end
        for s in (:title, :author, :ed, :fields, :notes)
            val = get(finfo, s, "")
            val != "" && setfield!(hdr, s, val)
        end
    end
    return 1
end

standardkind(s::AbstractString) = lowercase(replace(s, '-' => ' '))

function addsvd!(db::MatrixDatabase)
    added = 0
    for data in values(db.data)
        if data isa RemoteMatrixData && !issvdok(data)
            added += addsvd!(data)
        end
    end
    added
end

issvdok(data::RemoteMatrixData) = data.svdok
issvdok(::MatrixData) = false

"""
    downloadfile(url, out)

Copy file from remote or local url. Works around julia Downloads #69 and #36
"""
function downloadfile(url::AbstractString, out::AbstractString)
    wait(run(downloadcommand(url, out)))
    nothing
end
