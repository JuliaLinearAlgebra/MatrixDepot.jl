# return a list of file names without suffix in the directory
# e.g. filenames(mm) and filenames(uf)
function filenames(directory::String)

    namevec = String[]

    matdatadir = joinpath(Pkg.dir("MatrixDepot"), "data", "$directory")
    matvec = readdir(matdatadir)
    for file in matvec
        filename = split(file, '.')[1]
        push!(namevec, filename)
    end
    return namevec
end

# return a list of matrix data name in the collection
function matrix_data_name_list()
    matrices = String[]
    if isdir(joinpath(Pkg.dir("MatrixDepot"), "data", "uf"))
        for col in filenames("uf")
            for mat in filenames("uf/$(col)")
                push!(matrices, string(col, '/', mat))
            end
        end
    end

    if isdir(joinpath(Pkg.dir("MatrixDepot"), "data", "mm"))
        for col in filenames("mm")
            for d in filenames("mm/$(col)")
                for mat in filenames("mm/$(col)/$(d)")
                    push!(matrices, string(col, '/', d, '/', mat))
                end
            end
        end
    end
    matrices
end


# return a list of matrix name in the collection 
function matrix_name_list()
    matrices = sort(collect(keys(matrixdict)))
    append!(matrices, matrix_data_name_list())
    matrices
end

# return a list of groups in the collection
function group_list()
    groups = collect(keys(matrixclass))
    push!(groups, "data")
    groups = sort(groups)
    append!(groups, collect(keys(usermatrixclass)))
    groups
end

# print info about all matrices in the collection
function matrixdepot()
    # Print information strings
    println()
    println("Matrices:")

    matrices = matrix_name_list()

    i = 1
    for (index, mat) in enumerate(matrices)
        if i < 4 && length(mat) < 13
            i += 1
            @printf "%4d) %-14s" index mat
        else
            i = 1
            @printf "%4d) %-14s\n" index mat
        end
    end
    println()

    println("Groups:")

    groups = group_list()

    j = 1    
    for name in groups        
        if j < 4
            j += 1
            @printf "%12s" name
        else
            j = 1
            @printf "%12s\n" name
        end
    end
    println()
end

function matrixdepot{T}(name::String, ::Type{T}, m::Int, n::Int)
    # name is matrix name
    # m is the number of rows
    # n is the number of columns
    return matrixdict[name](T, m, n)
end
matrixdepot(name::String, m::Int, n::Int) = matrixdepot(name, Float64, m, n)

function matrixdepot{T}(name::String, ::Type{T}, n::Int, alpha)
    # name: matrix name
    # n: dimension of the matrix
    # alpha : scalar
    return matrixdict[name](T, n, alpha)
end
matrixdepot(name::String, n::Int, alpha) = matrixdepot(name, typeof(alpha), n, alpha)

function matrixdepot{T}(name::String, ::Type{T}, n::Int)
    # name is the matrix name
    # n is the dimension of the matrix (square)
    return matrixdict[name](T, n)
end

matrixdepot(name::String, n::Int) = matrixdepot(name, Float64, n)

function matrixdepot{T}(name::String, ::Type{T}, n::Int, alpha, beta)
    # name is the matrix name
    # n is the dimension of the matrix
    # alpha, beta are scalars
    return matrixdict[name](T, n, alpha, beta)
end
matrixdepot(name::String, n::Int, alpha, beta) = matrixdepot(name, Float64, n, alpha, beta)

function matrixdepot{T}(name::String, ::Type{T}, m::Int, n::Int, alpha, k::Int)
    # name is the matrix name
    # m, n are row and column dimensions of the matrix
    # alpha is a scalar
    # k is int
    return matrixdict[name](T, m, n, alpha, k)
end
matrixdepot(name::String, m::Int, n::Int, alpha, k::Int) = matrixdepot(name, Float64, m, n, alpha, k)

function matrixdepot{T}(name::String, ::Type{T}, m::Int, n::Int, alpha, theta)
    # name: matrix name
    # m, n are row and column dimensions of the matrix
    # alpha and theta are scalars
    return matrixdepot[name](T, m, n, alpha, theta)
end

function matrixdepot{T}(name::String, x::Vector{T}, y::Vector{T}, z::Vector{T})
    # name: matrix name
    # x, y, z: vectors
    return matrixdict[name](x,y,z)
end

function matrixdepot{T}(name::String, x::Vector{T}, y::Vector{T})
    # name: matrix name
    # x,y : vectors
    return matrixdict[name](x,y)
end


function matrixdepot{T}(name::String, x::Vector{T})
    return matrixdict[name](x)
end

function matrixdepot{T}(name::String, x::Vector{T}, n::Int)
    # name: matrix name
    # x: a vector
    # the column dimension of the matrix
    return matrixdict[name](x,n)
end


# Return information strings if name is a matrix name.
# Retuen a list of matrix names if name is a group.
function matrixdepot(name::String)
    # name is the matrix name or matrix properties
    if name in keys(matrixinfo)
        println(matrixinfo[name])
    elseif name in keys(matrixclass)
        return matrixclass[name]
    elseif name in keys(usermatrixclass)
        return usermatrixclass[name]
    elseif '/' in name  # print matrix data info
        namelist = split(name, '/')
        length(namelist) == 2 ? matdatadir = joinpath(Pkg.dir("MatrixDepot"), "data", "uf") :
                                matdatadir = joinpath(Pkg.dir("MatrixDepot"), "data", "mm")

        pathfilename = string(matdatadir, '/', name, ".mtx")
        println(ufinfo(pathfilename))
        println("use matrixdepot(\"$name\", :read) to read the data")
        return

    elseif name == "data" # deal with the group "data"
        return matrix_data_name_list()
    else
        error("\"$(name)\" is not included in Matrix Depot.")
    end
end

# access matrices by number
function matrixdepot(num::Int)
    matrixstrings = matrix_name_list()
    n = length(matrixstrings)
    if num > n
        error("There are $(n) parameterized matrices, but you ask for the $(num)-th ")
    end
    return matrixstrings[num]
end

function matrixdepot(I::UnitRange{Int})
    matrixnamelist = ASCIIString[]
    for i in I
        push!(matrixnamelist, matrixdepot(i))
    end
    return matrixnamelist
end

# generate the required matrix
# method = :read   (or :r) read matrix data
#          :get    (or :g) download matrix data
#          :search (or :s) search collection information
function matrixdepot(name::String, method::Symbol)
    if method == :r || method == :read
        length(split(name, '/')) == 2 ? matdatadir = joinpath(Pkg.dir("MatrixDepot"), "data", "uf") :
                                        matdatadir = joinpath(Pkg.dir("MatrixDepot"), "data", "mm")
        pathfilename = string(matdatadir, '/', name, ".mtx")

        if VERSION < v"0.4.0-dev+1419"
            return MatrixMarket.mmread(pathfilename)
        else
            return sparse(Base.SparseMatrix.CHOLMOD.Sparse(pathfilename))
        end
    elseif method == :g || method == :get
        MatrixDepot.get(name)
    elseif method == :s || method == :search
        MatrixDepot.search(name)
    else
        error("unknown symbol $method.
              use :read (or :r) to read matrix data;
              use :get  (or :g) to download matrix data;
              use :search (or :s) to search for collection information.")
    end
end


# Return a list of matrix names with common properties
# when multiple properties are given.
function matrixdepot(prop1::String, otherprops::String...)
    if length(otherprops) == 0
        matrixdepot(prop1)
    else
        commonprop = matrixdepot(prop1)
        for prop in otherprops
            commonprop = intersect(commonprop, matrixdepot(prop))
        end
    end
    return commonprop
end

#add new group
function addgroup(ex)
    propname = string(ex.args[1])
    !(propname in group_list()) || throw(ArgumentError("$propname is an existing group."))

    for matname in eval(ex.args[2])
        matname in matrix_name_list() || throw(ArgumentError("$matname is not in the collection."))
    end

    user = joinpath(Pkg.dir("MatrixDepot"), "src", "user.jl")
    s = readall(user)
    iofile = open(user, "w")
    newprop = s[1:end-4] * "\""  * propname * "\" => ["
    for str in eval(ex.args[2])
        newprop *= "\"" * str * "\", "
    end
    newprop = newprop * "],\n" * s[end-3:end]
    try
        write(iofile, newprop);
    finally
        close(iofile)
    end
end

macro addgroup(ex)
    esc(addgroup(ex))
end

# remove group
function rmgroup(ex)
    propname = string(ex)
    !(propname in keys(matrixclass)) || throw(ArgumentError("$propname can not be removed."))
    propname in keys(usermatrixclass) || throw(ArgumentError("Can not find group $propname."))

    user = joinpath(Pkg.dir("MatrixDepot"), "src", "user.jl")
    s = readall(user)
    iofile = open(user, "w")
    rg = Regex("""\"""" * eval(propname) * ".+")
    key = search(s, rg) # locate the propname in user.jl to remove.
    start_char = key[1] # the start of the line
    end_char = key[end] # the end of the line
    s = s[1:start_char - 1] * s[end_char+1:end]
    try
        write(iofile, s);
    finally
        close(iofile)
    end
end

macro rmgroup(ex)
    esc(rmgroup(ex))
end
