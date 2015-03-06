# return a list of file names in matrix database
function filenamevec()
    namevec = {}
    matdatadir = joinpath(Pkg.dir("MatrixDepot"), "data", "mat")
    matvec = readdir(matdatadir)
    for file in matvec
        filename = split(file, '.')[1]
        push!(namevec, filename)
    end
    return namevec
end

# print info about all matrices in the collection
function matrixdepot()
    # Print information strings 
    println()
    println("          | symmetric |  inverse  | ill-cond  |  pos-def  |   eigen   |")
    for mat in keys(matrixdict)
        @printf "%10s|" mat
        for prop in ["symmetric", "inverse", "ill-cond", "pos-def", "eigen"]
            if mat in matrixclass[prop]
                print("     *     |")
            else
                print("           |")
            end
        end
        println()
    end

    # Print UF sparse matrix files
    println()
    for file in filenamevec()
        @printf "%10s|" file
        print("  UF sparse matrix")
        println()
    end
    
    # print user defined properties
    if length(usermatrixclass) != 0
        println()
        println("New Properties:")
        println()
        for prop in keys(usermatrixclass)
            print(prop, " = [ ")
            for mat in usermatrixclass[prop]
                print(mat, ", ")
            end
            print("] \n")
            println()
        end
    end
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

function matrixdepot(name::String, n::Int)
    # name is the matrix name
    # n is the dimesion of the matrix
    # magic square, Pascal matrix and binomial matrix 
    # are exceptions: Int Array by default.
    if name in  ["magic", "pascal", "binomial"]
        matrixdepot(name, Int, n)
    else
        matrixdepot(name, Float64, n)
    end
end

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
# Retuen a list of matrix names if name is a property. 
function matrixdepot(name::String)
    # name is the matrix name or matrix properties
    if name in keys(matrixinfo)
        println(matrixinfo[name])
    elseif name in keys(matrixclass)
        return matrixclass[name]
    elseif name in keys(usermatrixclass)
        return usermatrixclass[name]
    elseif name in filenamevec()
        filename = string(name, ".mat")
        matdatadir = joinpath(Pkg.dir("MatrixDepot"), "data", "mat")
        pathfilename = string(matdatadir, '/', filename)
        matdata = matread(pathfilename)
        return (matdata["Problem"])["A"]
    else
        error("Your matrix or class is not included in Matrix Depot.")
    end
end

# return info for UF sparse matrices
function matrixdepot(name::String, s::Symbol)
    if s == :info
        filename = string(name, ".mat")
        matdatadir = joinpath(Pkg.dir("MatrixDepot"), "data", "mat")
        pathfilename = string(matdatadir, '/', filename)
        matdata = matread(pathfilename)
        return matdata["Problem"]
    else
        error("use symbol :info")
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

#addproperty
function addproperty(ex)
    propname = string(ex.args[1])
    !(propname in keys(matrixclass)) || throw(ArgumentError("$propname is an existing property."))
    !(propname in keys(usermatrixclass)) || throw (ArgumentError("You have defined property $propname."))
    for matname in eval(ex.args[2])
        matname in keys(matrixdict) || throw(ArgumentError("$matname is not in the collection."))
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

macro addproperty(ex)
    esc(addproperty(ex))
end

function rmproperty(ex)
    propname = string(ex)
    !(propname in keys(matrixclass)) || throw(ArgumentError("$propname can not be removed."))
    propname in keys(usermatrixclass) || throw (ArgumentError("Can not find property $propname."))
   
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

macro rmproperty(ex)
    esc(rmproperty(ex))
end
