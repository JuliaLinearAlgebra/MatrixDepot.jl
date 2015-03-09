# Download data from UF Sparse Matrix Collection and NIST Matrix Market

const UF_URL = "http://www.cise.ufl.edu/research/sparse/"  # UF Sparse Matrix collection
const MM_URL = "http://math.nist.gov/MatrixMarket/" # Matrix Market collectio
const DATA_DIR = joinpath(Pkg.dir("MatrixDepot"), "data")

# download html files and store matrix data as a list of tuples
function downloaddata(; collection::Symbol = :UF, generate_list::Bool = true)
    if collection == :UF     # UF Sparse matrix collection
        dlurl = string(UF_URL, "matrices/list_by_id.html")
        matrices = string(DATA_DIR, "/uf_matrices.html")
    elseif collection == :MM  # Matrix Market
        dlurl = string(MM_URL, "matrices.html")
        matrices = string(DATA_DIR, "/mm_matrices.html")
    else 
        error("unknown collection $collection")
    end
    isfile(matrices) || download(dlurl, matrices)
    if VERSION < v"0.4.0-dev+2197"
        matrixdata = {}
    else
        matrixdata = []
    end
   
    if generate_list
        open(matrices) do f
            if collection == :UF        
                for line in readlines(f)
                    
                    if contains(line, """MAT</a>""")
                        collectionname, matrixname = split(split(line, '"')[2], '/')[end-1:end]
                        matrixname = split(matrixname, '.')[1]
                        push!(matrixdata, (collectionname, matrixname)) 
                    end

                end
                
            elseif collection == :MM
                for line in readlines(f)

                    if contains(line, """<A HREF=\"/MatrixMarket/data/""")
                        collectionname, setname, matrixname = split(split(line, '"')[2], '/')[4:6]
                        matrixname = split(matrixname, '.')[1]
                        push!(matrixdata, (collectionname, setname, matrixname) )
                    end

                end 
            
            end
        
        end
    return matrixdata    
    end
   
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
    downloaddata(collection =:UF, generate_list = false)
    downloaddata(collection =:MM, generate_list = false)
end

# get
# --------------
# get(NAME) download a matrix from UF sparse matrix collection
# where NAME is a string of collection name + '/' + matrix name.
#
# Example
# -------
# get("HB/1138_bus", collection = :UF)
# get("Pajek/GD98_a")
# get("")
#
function get(name; collection::Symbol = :UF)
    if collection == :UF

    elseif collection == :MM

    else
        error("unknown collection $(collection)")
    end

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

