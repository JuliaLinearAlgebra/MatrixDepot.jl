########################
# helper functions
########################
using Markdown

# return a list of file names without suffix in the directory
# e.g. filenames(mm) and filenames(uf)
function filenames(directory::AbstractString)

    namevec = AbstractString[]

    matdatadir = joinpath(dirname(@__FILE__),"..", "data", "$directory")
    matvec = readdir(matdatadir)
    for file in matvec
        filename = rsplit(file, '.', limit=2)[1]
        push!(namevec, filename)
    end
    return namevec
end

# return the path to a folder inside the data directory
data_dir(name::AbstractString) = joinpath(dirname(@__FILE__), "..", "data", name)
user_file(name::AbstractString) = string(MY_DEPOT_DIR, "/$(name)")

# return a list of matrix data name in the collection
function matrix_data_name_list()
    matrices = AbstractString[]
    if isdir(data_dir("uf"))
        for col in filenames("uf")
            for mat in filenames("uf/$(col)")
                file = data_dir(joinpath("uf", col, mat, mat * ".mtx"))
                if isfile(file)
                    name = join((col, mat), "/")
                    push!(matrices, name)
                end
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
    matrices = sort(collect(keys(MATRIXDICT)))
    append!(matrices, matrix_data_name_list())
    matrices
end

# return a list of groups in the collection
_matrix_class() = collect(keys(MATRIXCLASS))
_user_matrix_class() = collect(keys(usermatrixclass))

function group_list()
    groups = _matrix_class()
    try
        append!(groups, _user_matrix_class())
    finally
        nothing
    end
    push!(groups, "data")
    push!(groups, "all")
    sort(groups)
end

##########################
# display information
##########################

# print info about all matrices in the collection
"""
`matrixdepot()`

Print all the matrices and groups in the collection.
"""
matrixdepot() = matrixdepot(stdout)

function matrixdepot(io::IO)
    # Print information strings
    println(io)
    println(io, "Matrices:")

    matrices = matrix_name_list()

    i = 1
    for (index, mat) in enumerate(matrices)
        if i < 4 && length(mat) < 14
            i += 1
            @printf io "%4d) %-14s" index mat
        else
            i = 1
            @printf io "%4d) %-14s\n" index mat
        end
    end
    println(io)

    println(io, "Groups:")

    groups = group_list()

    j = 1
    for name in groups
        if j < 4 && length(name) < 12
            j += 1
            @printf io "  %-12s" name
        else
            j = 1
            @printf io "  %-12s\n" name
        end
    end
    println(io)
end

#############################
# matrix generators
#############################

"""
`matrixdepot(matrix name, p1, p2...)`

Return a matrix specified by the query string `matrix name`.
`p1, p2...` are input parameters depending on `matrix name`.
"""
matrixdepot(name::AbstractString, args...) = MATRIXDICT[name](args...)

# generate the required matrix
# method = :read   (or :r) read matrix data
"""
`matrixdepot(data, symbol)`

Generate the data if `symbol = :r (or :read)`; download the data if `symbol = :g (or :get)`.
"""
function matrixdepot(name::AbstractString, method::Symbol, db::MatrixDatabase=MATRIX_DB;
                     meta::Bool = false)

    if method == :r || method == :read
        namelist = split(name, '/')
        if  length(namelist) == 2
            ufreader(string(data_dir("uf"), '/', namelist[1]), namelist[2], info = false, meta = meta)
        elseif length(namelist) == 3
            mmreader(data_dir("mm"), name, info = false)
        else
            matrixdepot(db.aliases[name], method, meta=meta)
        end
    else
        error("unknown symbol $method. use :read (or :r) to read matrix data")
    end
end


#########################
# access matrices
#########################


"""
`matrixdepot(number, range...)`

Access matrices by number, range or a mixture of numbers and ranges.
"""
function matrixdepot(num::Integer)
    matrixstrings = matrix_name_list()
    n = length(matrixstrings)
    if num > n
        throw(ArgumentError("There are $(n) parameterized matrices, but you asked for the $(num)-th."))
    end
    return matrixstrings[num]
end

function matrixdepot(ur::UnitRange)
    matrixnamelist = AbstractString[]
    for i in ur
        push!(matrixnamelist, matrixdepot(i))
    end
    return matrixnamelist
end

IntOrUnitRange = Union{Integer, UnitRange}

function matrixdepot(vs::IntOrUnitRange...)
    matrixnames = AbstractString[]
    for i in vs
        if typeof(i) <: Integer
            push!(matrixnames, matrixdepot(i))
        else
            append!(matrixnames, matrixdepot(i))
        end
    end
    return matrixnames
end

# Return a list of matrix names with common properties
# when multiple properties are given.
function matrixdepot(props::AbstractString...)
    common = matrixdepot(props[1])
    for prop in props[2:end]
        common = intersect(common, matrixdepot(prop))
    end
    return common
end


#######################
# matrix group
#######################

#add new group
function addgroup(ex)
    isdir(MY_DEPOT_DIR) || error("can not find directory '$MY_DEPOT_DIR'")
    propname = string(ex.args[1])
    !(propname in group_list()) || throw(ArgumentError("$propname is an existing group."))

    for matname in eval(ex.args[2])
        matname in matrix_name_list() || throw(ArgumentError("$matname is not in the collection."))
    end

    user = joinpath(user_file("group.jl"))
    s = read(user, String)
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
    !(propname in keys(MATRIXCLASS)) || throw(ArgumentError("$propname can not be removed."))
    propname in keys(usermatrixclass) || throw(ArgumentError("Can not find group $propname."))

    user = joinpath(user_file("group.jl"))
    s = read(user, String)
    iofile = open(user, "w")
    rg = Regex("""\"""" * eval(propname) * ".+")
    key = coalesce(findfirst(rg, s), 0:-1) # locate the propname in user.jl to remove.
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

abstract type MatrixGenerator end
abstract type FunctionName <: MatrixGenerator end
abstract type Group <: MatrixGenerator end


include_generator(::Type{FunctionName}, fn::AbstractString, f::Function) = (MATRIXDICT[fn] = f)

function include_generator(::Type{Group}, groupname::AbstractString, f::Function)
    if groupname in keys(MATRIXCLASS)
        push!(MATRIXCLASS[groupname], fname(f))
    elseif groupname in keys(usermatrixclass)
        push!(usermatrixclass[groupname], fname(f))
    else
        error("$(groupname) is not a group in MatrixDepot, use
              @addgroup to add this group")
    end
end

"return the name of the function `f` as a string."
function fname(f::Function)
    for (key, value) in MATRIXDICT
        if value == f
            return key
        end
    end
end

"return an array of full matrix names and aliases matching given pattern."
function list(r::Regex, db::MatrixDatabase=MATRIX_DB)
    result = AbstractString[]
    
    if startswith(r.pattern, "^#")
        for (alias, name) in db.aliases
            if match(r, alias) !== nothing
                push!(result, name)
            end
        end
    else
        for name in matrix_name_list()
            if match(r, name) !== nothing
                push!(result, name)
            end
        end
        for (alias, name) in db.aliases
            if !startswith(alias, '#') && match(r, alias) !== nothing && ! ( name in result )
                push!(result, name)
            end
            if match(r, name) !== nothing && ! ( name in result )
                push!(result, name)
            end
        end
    end
    sort!(result)
    unique!(result)
    result
end

function listdir(r::Regex, depth::Int, db::MatrixDatabase=MATRIX_DB)
    result = AbstractString[]
    f(x, n) = string(join(x[1:n], '/'), "/*" ^ max(length(x) - n, 0))
    for name in (f(x, depth) for x in split.(keys(db.data), '/') if length(x) >= depth)
        if match(r, name) !== nothing
            push!(result, name)
        end
    end
    sort!(result)
    unique!(result)
end

function list(p::AbstractString, db::MatrixDatabase=MATRIX_DB)
    p in keys(MATRIXCLASS) && ( return sort(MATRIXCLASS[p]) )
    p in keys(usermatrixclass) && ( return sort(usermatrixclass[p]) )
    p == "data" && ( return matrix_data_name_list() )
    p == "loaded" && ( return matrix_name_list() )
    p == "generated" && ( return sort!(collect(keys(MATRIXDICT))) )
    p == "all" && ( return sort!(union(matrix_name_list(), values(db.aliases))) )
   
    p = replace(p, "//+" => '/')
    depth = count(x == '/' for x in p)

    if occursin(r"[*?.]", p)
        p = replace(p, "**" => "([^/]+/)*[^/]+")
        p = replace(p, '*' => "[^/]+")
        p = replace(p, '?' => "[^/]")
        p = replace(p, '.' => "[.]")
    end
    if endswith(p, '/')
        length(p) == 1 && ( p = ".*/" )
        listdir(Regex(string('^', p)), depth, db)
    else
        
        list(Regex(string('^', p, '$')), db)
    end
end

flatten(a) = collect(Iterators.flatten(a))
list(i::Integer, db::MatrixDatabase=MATRIX_DB) = list(Regex(string('^', aliasid(i), '$')), db)
function list(r::OrdinalRange, db::MatrixDatabase=MATRIX_DB)
    listdb(r) = list(r, db)
    flatten(listdb.(aliasid.(r) ∩ keys(db.aliases)))
end
function list(r::AbstractVector, db::MatrixDatabase=MATRIX_DB)
    listdb(r) = list(r, db)
    flatten(listdb.(r))
end

const Pattern = Union{AbstractString,Regex,OrdinalRange,Integer,AbstractVector}

"return information about all matrices selected by pattern"
function info(p::Pattern, db::MatrixDatabase=MATRIX_DB)
    mdbuffer = Markdown.MD([])
    md = mdbuffer
    for name in list(p)
        try
            md = info(db.data[name])
            #=
            if data !== nothing && !isempty(data)
                display(data)
            end
            =#
        catch
            md = Markdown.parse("# $name\nno info available")
        finally
            append!(mdbuffer.content, md.content)
        end
    end
    mdbuffer
end

_mdheader(md::Markdown.MD, p, o) = isempty(md.content) ? (nothing, o) : _mdheader(md.content[1], md, o)
_mdheader(md::Markdown.Header, p, o) = (md, p)
_mdheader(md, p, o) = (nothing, o)

function info(data::GeneratedMatrixData)
    # func = data.fuc
    # md = @eval :(Docs.@doc $$func)
    md = eval(Meta.parse("Docs.@doc $(data.func)", raise = false))
    # As md is cached internally, need to make copies
    mdh, md = _mdheader(md, nothing, md)
    if mdh != nothing
        mdh, md = Markdown.Header(copy(mdh.text)), copy(md)
        push!(mdh.text, " ($(data.name))")
        md.content[1] = mdh
    end
    md
end

_repl(a::AbstractString) = a
_repl(a::AbstractString, p::Pair, q::Pair...) = _repl(replace(a, p), q...)

function info(data::RemoteMatrixData)
    txt = ufinfo(matrixfile(data))
    txt = _repl(txt, r"^%-+$"m => "---", r"^%%" => "###### ", r"%+" => "* ")
    md = Markdown.parse(txt)
    insert!(md.content, 1, Markdown.Header{1}([data.name]))
    md
end

"""
    matrix(p, [db,] args...)
return a matrix given a name
* `p` is a unique problem identifier (`String`, `Regex`, `Integer`)
* `db` is optional `MatrixDatabase` defaulting to global `MATRIX_DB`
* `args` are optional arguments - only used for generated matrices
"""
matrix(p::Pattern, db::MatrixDatabase=MATRIX_DB, args...) = matrix(mdopen(p, db), args...)
matrix(p::Pattern, args...) = matrix(mdopen(p), args...)

"""
    rhs(p[, db]
return right hand side of problem or `nothing`.
see also @matrix.
"""
rhs(p::Pattern, db::MatrixDatabase=MATRIX_DB) = rhs(mdopen(p))

"""
    solution(p[, db])
return solution of problem corresponding to right hand side or `nothing`.
see also @matrix.
"""
solution(p::Pattern, db::MatrixDatabase=MATRIX_DB) = solution(mdopen(p))

"""
    load(pattern[, db])
load data from remote repository for all problems matching pattern"
"""
function load(p::Pattern, db::MatrixDatabase=MATRIX_DB)
    for name in list(p)
        try
            loadmatrix(db.data[name], db)
        catch ex
            @warn "could not load $name: $ex"
        end
    end
end

"""
    mdopen(pattern[ ,db])
return `MatrixData` object, which can be used with data access functions.
If the pattern has not a unique resolution, en error is thrown.
"""
function mdopen(p::Pattern, db::MatrixDatabase=MATRIX_DB)
    li = list(p)
    length(li) == 0 && error("no matrix according to $p found")
    length(li) > 1  && error("pattern not unique: $p -> $li")
    db.data[li[1]]
end    

