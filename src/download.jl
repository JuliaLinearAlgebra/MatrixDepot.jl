
#####################################################
# Download data from UF Sparse Matrix Collection
#####################################################
const UF_URL = "http://www.cise.ufl.edu/research/sparse/"
const DATA_DIR = joinpath(Pkg.dir("MatrixDepot"), "data")

# download html files and store matrix data as a list of tuples
function downloaddata(; generate_list::Bool = true)
    # UF Sparse matrix collection
    dlurl = string(UF_URL, "matrices/list_by_id.html")
    matrices = string(DATA_DIR, "/uf_matrices.html")

    isfile(matrices) || download(dlurl, matrices)

    if generate_list
        matrixdata = Tuple[]
        open(matrices) do f

            for line in readlines(f)

                if contains(line, """MAT</a>""")
                    collectionname, matrixname = split(split(line, '"')[2], '/')[end-1:end]
                    matrixname = split(matrixname, '.')[1]
                    push!(matrixdata, (collectionname, matrixname))
                end

            end
        end
        return matrixdata
    end

end

# update database from the websites
function update()
    uf_matrices = string(DATA_DIR, "/uf_matrices.html")
    if isfile(uf_matrices)
        rm(uf_matrices)
    end

    downloaddata(generate_list = false)

end


function gunzip(fname)
    endswith(fname, ".gz") || error("gunzip: $fname: unknown suffix")

    destname = split(fname, ".gz")[1]

    open(destname, "w") do f
        GZip.open(fname) do g
            write(f, readall(g))
        end
    end
    destname
end


# get
# --------------
# get(NAME) download a matrix from UF sparse matrix collection
# where NAME is a string of collection name + '/' + matrix name.
#
# Example
# -------
# MatrixDepot.get("HB/1138_bus")
#
function get(name::String)

    if !isdir(string(DATA_DIR, '/', "uf"))
        mkdir(string(DATA_DIR, '/', "uf"))
    end
    matrixdata = downloaddata()
    collectionname, matrixname = split(name, '/')
    (collectionname, matrixname) in matrixdata ||
       error("can not find $collectionname/$matrixname in UF sparse matrix collection")
    fn = string(matrixname, ".tar.gz")
    uzfn = string(matrixname, ".mtx")
    url = string(UF_URL, "MM", '/', collectionname, '/', matrixname, ".tar.gz")

    dir = string(DATA_DIR, '/', "uf", '/', collectionname, '/')
    if !isdir(dir)
        mkdir(string(DATA_DIR, '/', "uf", '/', collectionname))
    end
    dirfn = string(dir, fn)
    diruzfn = string(dir, matrixname, '/', uzfn)


    !isfile(string(dir, uzfn)) || error("file $(uzfn) exits, no need to download")
    try
        download(url, dirfn)
        println("download:", dirfn)
    catch
        error("fail to download $fn")
    end


    if !isfile(diruzfn)
        gunzip(dirfn)
    end
    rm(dirfn)

    run(`tar -vxf $(dir)/$(matrixname).tar -C $(dir)`)
    cp("$(diruzfn)", "$(dir)/$(uzfn)")
    rm(string(dir,'/', matrixname, ".tar"))
    rm(string(dir, '/', matrixname), recursive=true)

end


