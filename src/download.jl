# Download data from the University of Florida Sparse Matrix Collection

const topurl = "http://www.cise.ufl.edu/research/sparse/"
const datadir = joinpath(Pkg.dir("MatrixDepot"), "data")

# download id.html and store matrixdata as a list of tuples
function downloaddata(dlurl = string(topurl, "matrices/list_by_id.html"))    
    matrices = string(datadir, "/matrices.html")
    isfile(matrices) || download(dlurl, matrices)
    matrixdata = {}
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

# downloadsparse
# --------------
# downloadsparse(NAME) download a matrix from UF sparse matrix collection
# where NAME is a string of collection name + '/' + matrix name.
#
# Example
# -------
# downloadsparse("HB/1138_bus")
# downloadsparse("Pajek/GD98_a")
#
function downloadsparse(name)
    matrixdata = downloaddata()
    collectionname, matrixname = split(name, '/')
    if (collectionname, matrixname) in matrixdata
        fn = string(matrixname, ".mat") # file name 
        dlfname = string(datadir, '/', "mat", '/', fn)
        if !isfile(dlfname)
            url = string(topurl, "mat", '/', collectionname, '/', "$fn")
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

