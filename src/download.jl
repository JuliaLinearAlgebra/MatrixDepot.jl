# The implementations are inspired by MatrixMarket.jl
# https://github.com/JuliaSparse/MatrixMarket.jl
# The MatrixMarket.jl package is licensed under the MIT Expat License:
# Copyright (c) 2013: Viral B. Shah.

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


#####################################################
# Download data from UF Sparse Matrix Collection
#####################################################
const UF_URL = "http://www.cise.ufl.edu/research/sparse/"
const MM_URL = "http://math.nist.gov/MatrixMarket/matrices.html"
const DATA_DIR = joinpath(dirname(@__FILE__),"..", "data")

# download html files and store matrix data as a list of tuples
function downloaddata(; generate_list::Bool = true)
    # UF Sparse matrix collection
    dlurl = string(UF_URL, "matrices/list_by_id.html")
    matrices = string(DATA_DIR, "/uf_matrices.html")
    mm_matrices = string(DATA_DIR, "/mm_matrices.html")

    isfile(matrices) || download(dlurl, matrices)
    isfile(mm_matrices) || download(MM_URL, mm_matrices)

    if generate_list
        uf_matrixdata = Tuple[]
        mm_matrixdata = Tuple[]
        open(matrices) do f
            for line in readlines(f)
                if contains(line, """MAT</a>""")
                    collectionname, matrixname = split(split(line, '"')[2], '/')[end-1:end]
                    matrixname = split(matrixname, '.')[1]
                    push!(uf_matrixdata, (collectionname, matrixname))
                end
            end
        end
        
        open(mm_matrices) do f
            for line in readlines(f)
                if contains(line, """<A HREF="/MatrixMarket/data/""")
                    collectionname, setname, matrixname = split(split(line, '"')[2], '/')[4:6]
                    matrixname = split(matrixname, '.')[1]
                    push!(mm_matrixdata, (collectionname, setname, matrixname))
                end
            end
        end
        return uf_matrixdata, mm_matrixdata
    end

end

# given a matrix name, 
function search(matrixname::AbstractString)
    uf_matrixdata, mm_matrixdata = downloaddata()
    uf_matrices = AbstractString[]
    mm_matrices = AbstractString[]
    [push!(uf_matrices, m[2]) for m in uf_matrixdata]
    [push!(mm_matrices, m[3]) for m in mm_matrixdata]

    datalist = AbstractString[]
    if (matrixname in uf_matrices) || (matrixname in mm_matrices)
        uf_index = findin(uf_matrices, [matrixname])
        mm_index = findin(mm_matrices, [matrixname])
        if uf_index != 0
            for i in uf_index
                collectionname, matrixname = uf_matrixdata[i]
                push!(datalist, string(collectionname, '/', matrixname))
            end
        end
        if mm_index != 0
            for i in mm_index
                collectionname, setname, matrixname = mm_matrixdata[i]
                push!(datalist, string(collectionname, '/', setname, '/', matrixname))
            end
        end
        return datalist
    else
        error("can not find $matrixname in the database.")
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
# MatrixDepot.get("HB/1138_bus") # uf sparse matrix
# MatrixDepot.get("1138_bus") # Matrix Market
#
function get(name::AbstractString)   

    isdir(string(DATA_DIR,'/', "uf")) || mkdir(string(DATA_DIR, '/', "uf"))
    isdir(string(DATA_DIR,'/', "mm")) || mkdir(string(DATA_DIR, '/', "mm"))
                                               
    uf_matrixdata, mm_matrixdata = downloaddata()

    namelist = split(name, '/')
    collectionname = namelist[1]
    if length(namelist) == 2 # UF sparse matrix collection
        matrixnames = AbstractString[]
        if namelist[2] == "*"
            println("Downloading all matrices in group $(collectionname)...")
            for uf in uf_matrixdata
                if uf[1] == collectionname
                    push!(matrixnames, uf[2])
                end
            end
        else
            push!(matrixnames, namelist[2])
        end
        for matrixname in matrixnames           
            (collectionname, matrixname) in uf_matrixdata ||
                       error("can not find $collectionname/$matrixname in UF sparse matrix collection")
            fn = string(matrixname, ".tar.gz")
            uzfn = string(matrixname, ".mtx")
            url = string(UF_URL, "MM", '/', collectionname, '/', matrixname, ".tar.gz")
    
            dir = string(DATA_DIR, '/', "uf", '/', collectionname)
            if !isdir(dir)
                mkdir(dir)
            end
            
            dirfn = string(dir, '/', fn)
            diruzfn = string(dir, '/', matrixname, '/', uzfn)
            
            if isfile(diruzfn)
                error("$(collectionname)/$(matrixname) exits, no need to download")
            end

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
            rm(string(dir,'/', matrixname, ".tar"))
        end
    elseif length(namelist) == 3 # Matrix Market collection 
        
        collectionname, setname, matrixname = namelist
        mtxfname = string(matrixname, ".mtx")
        uzfn = string(matrixname, ".mtx.gz")
        
        dir = string(DATA_DIR, '/', "mm", '/', collectionname)
        isdir(dir) || mkdir(string(DATA_DIR, '/', "mm", '/', collectionname))
   
        dir = string(DATA_DIR, '/', "mm", '/', collectionname, '/', setname)
        isdir(dir) || mkdir(string(DATA_DIR, '/', "mm", '/', collectionname, '/', setname))
        
        dirfn = string(dir, '/', mtxfname)
        diruzfn = string(dir, '/', uzfn)
        
        !isfile(dirfn) || error("file $(mtxfname) exits, no need to download")
        url =  "ftp://math.nist.gov/pub/MatrixMarket2/$(collectionname)/$(setname)/$(matrixname).mtx.gz"

        try 
            download(url, diruzfn)
            println("download:", diruzfn)
        catch
            error("fail to download $uzfn")
        end
        gunzip(diruzfn)
        rm(diruzfn)

    elseif length(namelist) == 1
        stringvec = search(name)
        if length(stringvec) == 1
            return matrixdepot(stringvec[1], :get)
        else
            println("Try MatrixDepot.get(`name`), where `name` is one of the elements in the following Array:")
            return stringvec
        end
    else
        error("can not find $name")
    end
end
