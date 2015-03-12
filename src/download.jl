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


function gunzip(fname, newname)
    endswith(fname, ".gz") || error("gunzip: $fname: unknown suffix")
    
    fn = string(newname, ".mtx")
    open(fn, "w") do f
        GZip.open(fname) do g
            write(f, readall(g))
        end
    end
    fn
end


# get
# --------------
# get(NAME) download a matrix from UF sparse matrix collection
# where NAME is a string of collection name + '/' + matrix name.
# 
# get(NAME, collection = :MM) download a matrix from Matrix Market
# where NAME is string of collection name + '/' + set name + '/' + matrix name.
#
# Example
# -------
# MatrixDepot.get("HB/1138_bus", collection = :UF)
# MatrixDepot.get("Pajek/GD98_a")
# MatrixDepot.get("SPARSKIT/fidap/fidap020", collection =:MM)
#
function get(name::String, newname::String; collection::Symbol = :UF)
         
    if collection == :UF
        matrixdata = downloaddata()
        collectionname, matrixname = split(name, '/')
        (collectionname, matrixname) in matrixdata || 
                            error("can not find $collectionname\$matrixname in UF sparse matrix collection")
        fn = string(newname, ".mat")
       
        url = string(UF_URL, "mat", '/', collectionname, '/', matrixname, ".mat")
        dirfn = string(DATA_DIR, '/', "uf",'/',fn)
        old_dirfn = dirfn
    elseif collection == :MM
        matrixdata = downloaddata(collection =:MM)
        collectionname, setname, matrixname = split(name, '/')
        (collectionname, setname, matrixname) in matrixdata ||
                            error("can not find $collectionname/$setname/$matrixname in Matrix Market")
        fn = string(newname, ".mtx.gz")
        url = "ftp://math.nist.gov/pub/MatrixMarket2/$collectionname/$setname/$matrixname.mtx.gz"
        dirfn = string(DATA_DIR, '/',"mm", '/', fn)
        old_dirfn = string(DATA_DIR, '/',"mm", '/', newname, ".mtx" )
    else
        error("unknown collection $(collection)")
    end

   
    !isfile(ord_dirfn) || error("file $fn exists, try rename it with MatrixDepot.get(name, newname, collection)")    
    try 
        download(url, dirfn)
    catch
        error("fail to download $fn")
    end
    
    
    if collection == :MM
        fn = string(newname, ".mtx")
        if !isfile(fn) 
            gunzip(dirfn, newname)
        end
        rm(dirfn)
    end
    
end

function get(name::String; collection::Symbol = :UF)
    if collection == :UF
        collectionname, matrixname = split(name, '/')
        get(name, matrixname)
    elseif collection == :MM
        collectionname, setname, matrixname = split(name, '/')
        get(name, matrixname, collection = :MM)
    end
end 

