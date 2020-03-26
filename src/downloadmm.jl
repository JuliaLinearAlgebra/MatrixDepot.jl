
#####################################################
# Download data from MatrixMarket Collection
#####################################################

function readindex(remote::RemoteType, db::MatrixDatabase)
    file = downloadindex(remote)
    extract_names(db, file)
end

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
    nothing
end

# read indexfile
function downloadindex(remote::RemoteType)
    file = localindex(remote)
    url = redirect(indexurl(remote))
    if !isfile(file)
        println("dowloading index file $url")
        download(url, file)
    end
    file
end

# name translations (most of the MM problems are found with similar name in SuiteSparse)
function namemm2ss(mmname::AbstractString)
    s = split(mmname, '/')
    length(s) == 3 || return mmname
    s1, s2, s3 = s
    if s2 == "qcd"
        replace(s3, r"^([^.]*)\.(.)-00l(......)00$" => s"QCD/\1_\2-\3")
    elseif s2 == "tokamak"
        string("TOKAMAK/", s3)
    elseif s2 == "fidap"
        replace(s3, r"fidap0*([^0]*)" => s"FIDAP/ex\1")
    elseif s1 == "Harwell-Boeing"
        string("*/", replace(s3, r"_+" => "_"))
    elseif s1 == "NEP"
        if s2 == "crystal"
            string("*/", replace(s3, r"^(cry)(.*)$" => s"\1g\2"))
        elseif s3 == "rdb2048l"
            string("*/", "rdb2048")
        elseif s3 == "rdb2048"
            string("*/", "rdb2048_noL")
        else
            string("*/", replace(s3, r"^(bfw|mhd|rbs|odep)([\d]{2,3})([a-z])$" => s"\1\3\2"))
        end
    else
        string("*/", s3)
    end
end

function namess2mm(ssname::AbstractString)
    s = split(ssname, '/')
    length(s) == 2 || return ssname
    s1, s2 = s
    if s1 == "QCD"
        replace(s2, r"^([^_]*)_(.)-(...)-(....)$" => s"misc/qcd/\1.\2-00l\3-\4") * "00"
    elseif s1 == "TOKAMAK"
        string("SPARSKIT/tokamak", s2)
    elseif s1 == "FIDAP"
        x = replace(s2, r"^ex([0-9]*)$" => s"\1")
        x == s2 ? string(s1, '/', x) : "SPARSKIT/fidap/fidap$(printfint(parse(Int,x), 3))"
    elseif s1 == "HB"
        n = length(s2)
        N = startswith(s2, "watt_") || startswith(s2, "pores_") ? 7 : 8
        y = if n >= N
            s2
        else
            x = replace(s2, r"^(.*)(_+)(.*[0-9])$" => s"\1/\3")
            x == s2 ? s2 : replace(x, r"/" => "_" ^ (N + 1 - n))
        end
        "Harwell-Boeing/*/" * y
    elseif s1 == "Bai" && startswith(s2, "cryg")
        replace(s2, r"^(cry)g(.*)$" => s"NEP/crystal/\1\2")
    elseif s1 == "Bai"
        x = (s2 == "rdb2048" ? "rdb2048l" : s2 == "rdb2048_noL" ? "rdb2048" : s2)
        x = replace(x, r"^(bfw|mhd|rbs|odep)([a-z])([\d]{2,3})" => s"\1\3\2")
        string("NEP/*/", x)
    else
        string("*/*/", s2)
    end
end

printfint(n::Integer, pad::Int) = String(reverse(digits(n, base=UInt8(10), pad=pad)) .+ '0')

