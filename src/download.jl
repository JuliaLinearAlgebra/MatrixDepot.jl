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
const MM_URL2 = "ftp://math.nist.gov/pub/MatrixMarket2"
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
                if occursin("""MAT</a>""", line)
                    collectionname, matrixname = split(split(line, '"')[2], '/')[end-1:end]
                    matrixname = split(matrixname, '.')[1]
                    push!(uf_matrixdata, (collectionname, matrixname))
                end
            end
        end
        
        open(mm_matrices) do f
            for line in readlines(f)
                if occursin("""<A HREF="/MatrixMarket/data/""", line)
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
        uf_index = findall(isequal(matrixname), uf_matrices)
        mm_index = findall(isequal(matrixname), mm_matrices)
        for i in uf_index
            collectionname, matrixname = uf_matrixdata[i]
            push!(datalist, string(collectionname, '/', matrixname))
        end
        for i in mm_index
            collectionname, setname, matrixname = mm_matrixdata[i]
            push!(datalist, string(collectionname, '/', setname, '/', matrixname))
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
            write(f, read(g))
        end
    end
    destname
end

function matchnames(p, a)
    n = length(a)
    length(p) == n || (return false)
    for i = 1:n
        p[i] in ("*", a[i]) || (return false)
    end
    true
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
    if length(namelist) == 2 # UF sparse matrix collection
        matrixnames = Tuple{<:AbstractString,<:AbstractString}[]
        println("Downloading all matrices named $name...")
        for uf in uf_matrixdata
            if matchnames(namelist, uf)
                push!(matrixnames, uf)
            end
        end
        for fullname in matrixnames           
            collectionname, matrixname = fullname
            fn = string(matrixname, ".tar.gz")
            uzfn = string(matrixname, ".mtx")
            url = join((UF_URL, "MM", fullname...), "/") * ".tar.gz"

            dir = joinpath(DATA_DIR, "uf", collectionname)
            isdir(dir) || mkdir(dir)
            
            dirfn = joinpath(dir, fn)
            diruzfn = joinpath(dir, matrixname, uzfn)
            
            if isfile(diruzfn)
                continue
            end

            try
                println("downloading: ", url)
                download(url, dirfn)
            catch
                error("failed to download $uzfn")
            end
            
            if !isfile(diruzfn)
                gunzip(dirfn)
            end
            rm(dirfn)
            run(`tar -vxf $(dir)/$(matrixname).tar -C $(dir)`)
            rm(string(dir,'/', matrixname, ".tar"))
        end
    elseif length(namelist) == 3 # Matrix Market collection 
        
        matrixnames = Tuple{<:AbstractString,<:AbstractString,<:AbstractString}[]
        println("Downloading all matrices named $name...")
        for mm in mm_matrixdata
            if matchnames( namelist, mm)
                push!(matrixnames, mm)
            end
        end
        
        for fullname in matrixnames
            collectionname, setname, matrixname = fullname
            mtxfname = string(matrixname, ".mtx")
            uzfn = string(matrixname, ".mtx.gz")
            
            dir = joinpath(DATA_DIR, "mm", collectionname, setname)
            isdir(dir) || mkpath(dir)
            
            dirfn = joinpath(dir, mtxfname)
            diruzfn = joinpath(dir, uzfn)
            
            !isfile(dirfn) || continue
            url = join((MM_URL2, fullname...), "/") * ".mtx.gz"

            try 
                println("downloading: ", url)
                download(url, diruzfn)
            catch
                error("failed to download $uzfn")
            end
            gunzip(diruzfn)
            rm(diruzfn)
        end

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
