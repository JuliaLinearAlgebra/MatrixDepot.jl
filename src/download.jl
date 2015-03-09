# Download data from the University of Florida Sparse Matrix Collection

const topurl = "http://www.cise.ufl.edu/research/sparse/"
const datadir = joinpath(Pkg.dir("MatrixDepot"), "data")

# download id.html and store matrixdata as a list of tuples
function downloadUFdata(dlurl = string(topurl, "matrices/list_by_id.html"))    
    matrices = string(datadir, "/matrices.html")
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
    matrices = string(datadir, "/matrices.html")
    if isfile(matrices)
        rm(matrices)
        download(string(topurl, "matrices/list_by_id.html"), matrices)
    end        
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

    # check if the matrix is in the database
    if (collectionname, matrixname) in matrixdata  
        fn = string(matrixname, ".mat")              # file name 
        dlfname = string(datadir, '/', "mat", '/', collectionname, '_' ,fn)
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

function downloadMM(name)

end
