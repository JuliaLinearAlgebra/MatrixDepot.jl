# return info comment strings for UF sparse matrix
function ufinfo(filename::AbstractString)
    mmfile = open(filename,"r")
    info = "\n"
    ll = readline(mmfile)
    while length(ll) > 0 && ll[1] == '%'
        info = string(info, ll)
        ll = readline(mmfile)
    end
    return info
end

# read Matrix Market data
"""
`mmreader(dir, name, info)`

dir: directory of the file
name: file name
info: whether to return infomation
"""
function mmreader(dir::AbstractString, name::AbstractString; info::Bool = true)
    pathfilename = string(dir, '/', name, ".mtx")
    if info
        println(ufinfo(pathfilename))
        println("use matrixdepot(\"$name\", :read) to read the data")
    else
        sparse(Base.SparseMatrix.CHOLMOD.Sparse(pathfilename))
    end
end
    
# read UF sparse matrix data
"""
`ufreader(dir, name, info, meta)`

dir: directory of the file
name: file name
info: whether to return information
meta: wehterh to return metadata
"""
function ufreader(dir::AbstractString, name::AbstractString;
                  info::Bool = true, meta::Bool = false)
    dirname = string(dir, '/', name)
    files = filenames(dirname)
    if info
        println(ufinfo(string(dirname, '/', name, ".mtx")))
        if length(files) > 1
            println("metadata $(files)")
        end
        #println("use matrixdepot(\"$name\", :read) to read the data")
    else
        A = sparse(Base.SparseMatrix.CHOLMOD.Sparse(string(dirname, '/', name, ".mtx")))
        if meta
            metadict = Dict{AbstractString, Any}()
            datafiles = readdir(dirname)
            for data in datafiles
                dataname = split(data, '.')[1]
                try
                    metadict[dataname] =  sparse(Base.SparseMatrix.CHOLMOD.Sparse(string(dirname, '/', data)))
                catch
                    metadict[dataname] = readall(string(dirname,'/', data))
                end
            end
            metadict
        else
            A
        end        
    end
end
