########################
# helper functions
########################

# return a list of file names without suffix in the directory
# e.g. filenames(mm) and filenames(uf)
function filenames(directory::AbstractString)

    namevec = AbstractString[]

    matdatadir = joinpath(dirname(@__FILE__),"..", "data", "$directory")
    matvec = readdir(matdatadir)
    for file in matvec
        filename = split(file, '.')[1]
        push!(namevec, filename)
    end
    return namevec
end

# return the path to a folder inside the data directory
data_dir(name::AbstractString) = joinpath(dirname(@__FILE__), "..", "data", name)
user_file(name::AbstractString) = joinpath(dirname(@__FILE__), "..", "user", name)

# return a list of matrix data name in the collection
function matrix_data_name_list()
    matrices = AbstractString[]
    if isdir(data_dir("uf"))
        for col in filenames("uf")
            for mat in filenames("uf/$(col)")
                push!(matrices, string(col, '/', mat))
            end
        end
    end

    if isdir(data_dir("mm"))
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
    push!(groups, "all")
    groups = sort(groups)
    append!(groups, collect(keys(usermatrixclass)))
    groups
end

##########################
# display information
##########################

# print info about all matrices in the collection
"""
`matrixdepot()` 

Print all the matrices and groups in the collection.
"""
function matrixdepot()
    # Print information strings
    println()
    println("Matrices:")

    matrices = matrix_name_list()

    i = 1
    for (index, mat) in enumerate(matrices)
        if i < 4 && length(mat) < 14
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
        if j < 4 && length(name) < 12
            j += 1
            @printf "  %-12s" name
        else
            j = 1
            @printf "  %-12s\n" name
        end
    end
    println()
end

# Return information strings if name is a matrix name
# and return a list of matrix names if name is a group.
"""
`matrixdepot(name)` 

Return the documentation if `name` is a matrix name;
return a list of matrix names if `name` is a group name.
"""
function matrixdepot(name::AbstractString)
    # name is the matrix name or matrix properties
    if name in keys(matrixinfo)
        println(matrixinfo[name])
    elseif name in keys(matrixclass)
        matrices = matrixclass[name]
        return sort(matrices)
    elseif name in keys(usermatrixclass)
        matrices = usermatrixclass[name]
        return sort(matrices)
    elseif '/' in name  # print matrix data info
        namelist = split(name, '/')
        length(namelist) == 2 ? matdatadir = data_dir("uf") :
                                matdatadir = data_dir("mm")

        pathfilename = string(matdatadir, '/', name, ".mtx")
        println(ufinfo(pathfilename))
        println("use matrixdepot(\"$name\", :read) to read the data")
        return

    elseif name == "data" # deal with the group "data"
        return matrix_data_name_list()
    elseif name == "all" # all the matrix names in the collection
        return matrix_name_list()
    else
        throw(ArgumentError("No information is available for \"$(name)\"."))
    end
end


#############################
# matrix generators 
#############################

"""
`matrixdepot(matrix name, p1, p2...)` 

Return a matrix specified by the query string `matrix name`. 
`p1, p2...` are input parameters depending on `matrix name`.
"""
matrixdepot{T}(name::AbstractString, ::Type{T}, args...) = length(args) == 1 ? 
             matrixdict[name](T, args[1]) : matrixdict[name](T, collect(args)...)
matrixdepot(name::AbstractString, args...) = matrixdepot(name::AbstractString, Float64, args...)

# generate the required matrix
# method = :read   (or :r) read matrix data
#          :get    (or :g) download matrix data
#          :search (or :s) search collection information
"""
`matrixdepot(data, symbol)` 

Generate the data if `symbol = :r (or :read)`; download the data if `symbol = :g (or :get)`.
"""
function matrixdepot(name::AbstractString, method::Symbol)
    if method == :r || method == :read
        length(split(name, '/')) == 2 ? matdatadir = data_dir("uf"):
                                        matdatadir = data_dir("mm")
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
        throw(ArgumentError("unknown symbol $method.
              use :read (or :r) to read matrix data;
              use :get  (or :g) to download matrix data;
              use :search (or :s) to search for collection information."))
    end
end


#########################
# access matrices
#########################

# access matrices by number
"""
`matrixdepot(number, range...)` 

Access matrices by number, range or a mixture of numbers and ranges.
"""
function matrixdepot(num::Int)
    matrixstrings = matrix_name_list()
    n = length(matrixstrings)
    if num > n
        error("There are $(n) parameterized matrices, but you ask for the $(num)-th ")
    end
    return matrixstrings[num]
end

function matrixdepot(I::UnitRange{Int})
    matrixnamelist = AbstractString[]
    for i in I
        push!(matrixnamelist, matrixdepot(i))
    end
    return matrixnamelist
end

IntOrUnitRange = @compat Union{Int, UnitRange{Int}}

function matrixdepot(vs::IntOrUnitRange...)
    matrixnames = AbstractString[]
    for i in vs
        if typeof(i) <: Int
            push!(matrixnames, matrixdepot(i))
        else
            append!(matrixnames, matrixdepot(i))
        end
    end
    return matrixnames
end

# Return a list of matrix names with common properties
# when multiple properties are given.
function matrixdepot(prop1::AbstractString, otherprops::AbstractString...)
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


#######################
# matrix group
#######################

#add new group
function addgroup(ex)
    propname = string(ex.args[1])
    !(propname in group_list()) || throw(ArgumentError("$propname is an existing group."))

    for matname in eval(ex.args[2])
        matname in matrix_name_list() || throw(ArgumentError("$matname is not in the collection."))
    end

    user = joinpath(user_file("group.jl"))
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

"add a group to Matrix Depot"
macro addgroup(ex)
    esc(addgroup(ex))
end

# remove an added group
function rmgroup(ex)
    propname = string(ex)
    !(propname in keys(matrixclass)) || throw(ArgumentError("$propname can not be removed."))
    propname in keys(usermatrixclass) || throw(ArgumentError("Can not find group $propname."))

    user = joinpath(user_file("group.jl"))
    s = readall(user)
    iofile = open(user, "w")
    rg = Regex("""\"""" * eval(propname) * ".+")
    key = search(s, rg) # locate the propname in user.jl to remove.
    start_char = key[1] # the start of the line
    end_char = key[end] # the end of the line
    s = s[1:start_char - 2] * s[end_char+1:end]
    try
        write(iofile, s);
    finally
        close(iofile)
    end
end

"remove an added group from Matrix Depot"
macro rmgroup(ex)
    esc(rmgroup(ex))
end

################################
# user defined matrix generators
################################

abstract MatrixGenerator

abstract FunctionName <: MatrixGenerator
abstract Help <: MatrixGenerator
abstract Group <: MatrixGenerator


include_generator(::Type{FunctionName}, fn::AbstractString, f::Function) = (matrixdict[fn] = f)
include_generator(::Type{Help}, helplines::AbstractString, f::Function) = (matrixinfo[fname(f)] = helplines)
function include_generator(::Type{Group}, groupname::AbstractString, f::Function) 
    if groupname in keys(matrixclass)
        push!(matrixclass[groupname], fname(f))
    elseif groupname in keys(usermatrixclass)
        push!(usermatrixclass[groupname], fname(f))
    else
        error("$(groupname) is not a group in MatrixDepot, use
              @addgroup to add this group")
    end
end

"return the name of the function `f` as a string."
fname(f::Function) = split(string(f), '.')[2]
