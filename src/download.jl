# Download data from the University of Florida Sparse Matrix Collection

const UF_URL = "http://www.cise.ufl.edu/research/sparse/"
const DATA_DIR = joinpath(Pkg.dir("MatrixDepot"), "data")

# download id.html and store matrixdata as a list of tuples
function downloaddata(dlurl = string(UF_URL, "matrices/list_by_id.html"))    
    matrices = string(DATA_DIR, "/matrices.html")
    isfile(matrices) || download(dlurl, matrices)
    if VERSION < v"0.4.0-dev+2197"
        matrixdata = {}
    else
        matrixdata = []
    end
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

# update database from the website
function updatesparse()
    matrices = string(DATA_DIR, "/matrices.html")
    if isfile(matrices)
        rm(matrices)
        download(string(UF_URL, "matrices/list_by_id.html"), matrices)
    end        
end

# downloadUF
# --------------
# downloadUF(NAME) download a matrix from UF sparse matrix collection
# where NAME is a string of collection name + '/' + matrix name.
#
# Example
# -------
# downloadUF("HB/1138_bus")
# downloadUF("Pajek/GD98_a")
#
function downloadsparse(name; collection::Symbol = :UF)
    matrixdata = downloaddata()
    collectionname, matrixname = split(name, '/')

    # check if the matrix is in the database
    if (collectionname, matrixname) in matrixdata  
        fn = string(matrixname, ".mat")              # file name 
        dlfname = string(DATA_DIR, '/', "mat", '/', collectionname, '_' ,fn)
        if !isfile(dlfname)
            url = string(UF_URL, "mat", '/', collectionname, '/', "$fn")
            try 
                download(url, dlfname)
            catch
                error("fail to download $fn")
            end
        end
    else
        error("can not find $collectionname/$matrixname in database")
    end
end

