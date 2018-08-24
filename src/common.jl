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

function group_list()
    groups = Symbol[]
    append!(groups, keys(MATRIXCLASS))
    append!(groups, keys(usermatrixclass))
    append!(groups, keys(SUBSETS))
    sort(groups)
end

##########################
# display information
##########################

# print info about all matrices in the collection
"""
    overview([db])

return formatted overview about matrices and groups in the collection.
"""
function overview(db::MatrixDatabase=MATRIX_DB)
    # Print information strings
    io = IOBuffer()
    println(io)
    println(io, "Matrices:")

    matrices = list(:local)

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
        if j < 4 && length(string(name)) < 12
            j += 1
            @printf io "  %-12s" name
        else
            j = 1
            @printf io "  %-12s\n" name
        end
    end
    println(io)
    String(take!(io))
end

#######################
# matrix group
#######################

# write one property association
function propline(io::IO, propname, matnames)
    println(matnames)
    write(io, repr(propname))
    write(io, " => [")
    for str in matnames
        write(io, repr(str))
        write(io, ", ")
    end
    write(io, "],\n")
end

# add, remove, or replace complete user group
function modgroup(prop::Symbol, mats::Union{Nothing,Vector{<:AbstractString}})
    prop in keys(MATRIXCLASS) && throw(ArgumentError("$prop can not be modified."))

    user = abspath(MY_DEPOT_DIR, "group.jl")
    s = read(user, String)          # read complete file into s
    rg = Regex(repr(prop) * r"\W*=>\W*(\[.*\]\W*,\W*\n)".pattern)
    ppos = findfirst(rg, s)         # locate the prop in user.jl to remove.
    if ppos !== nothing
        start_char = first(ppos) - 1    # the start of the line
        end_char = last(ppos)           # the end of the line
    else
        ppos = findnext(r"\);", s, 1)
        if ppos !== nothing
            start_char = end_char = first(ppos) - 1
        else
            start_char = end_char = length(s)
        end
    end
    if mats !== nothing
        mats = sort(mats)
    end
    open(user, "w") do io
        write(io, s[1:start_char])
        if mats !== nothing
            propline(io, prop, mats)
        end
        write(io, s[end_char+1:end])
    end
    if mats !== nothing
        usermatrixclass[prop] = mats
    else
        delete!(usermatrixclass, prop)
    end
    return nothing
end

"add a group to Matrix Depot"
macro addgroup(ex)
    name = Symbol(ex.args[1])
    esc(modgroup(name, eval(ex.args[2])))
end

"add or replace group in Matrix Depot"
macro modifygroup(ex)
    name = Symbol(ex.args[1])
    esc(modgroup(name, eval(ex.args[2])))
end

"remove an added group from Matrix Depot"
macro rmgroup(ex)
    name = Symbol(ex)
    esc(modgroup(name, nothing))
end

################################
# user defined matrix generators
################################

abstract type MatrixGenerator end
abstract type FunctionName <: MatrixGenerator end
abstract type Group <: MatrixGenerator end


function include_generator(::Type{FunctionName}, fn::AbstractString, f::Function)
    (haskey(MATRIXDICT, fn) ? MATRIXDICT : USERMATRIXDICT)[fn] = f
end

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
        value == f && return key
    end
    for (key, value) in USERMATRIXDICT
        value == f && return key
    end
    "unknown-function"
end

"""
    list(p::Pattern[, db]
return a vector of full matrix names where name or alias match given pattern.
`p` can be one of the following:
+ a plain string (without characters `*` and `?`) which must match exactly
+ a string containing `*` and `?` acting like a shell path pattern
+ a regular expression
+ an integer matching equivalent to the alias string `"#\$p"`
+ a range of integers
+ one of the symbols `:local`, `:remote`, `:loaded`, `:all`
+ a vector of patterns meaning the union
+ a tuple of patterns meaning the intersection
"""
function list(r::Regex, db::MatrixDatabase=MATRIX_DB)
    result = AbstractString[]
    
    if startswith(r.pattern, "^#")
        for (alias, name) in db.aliases
            if match(r, alias) !== nothing
                push!(result, name)
            end
        end
    else
        for name in keys(db.data)
            if match(r, name) !== nothing
                push!(result, name)
            end
        end
    end
    sort!(result)
    unique!(result)
    result
end

function listdir(r::Regex, depth::Int, db::MatrixDatabase=MATRIX_DB)
    result = Dict{AbstractString, Int}()
    f(x, n) = string(join(x[1:n], '/'), "/*" ^ max(length(x) - n, 0))
    for name in list(r)
        li = split(name, '/')
        if length(li) >= depth
            key = f(li, depth)
            result[key] = get!(result, key, 0) + 1
        end
    end
    sort!([string(k, " - (", v, ")") for (k,v) in result])
end

function list(p::Symbol, db::MatrixDatabase=MATRIX_DB)
    if haskey(SUBSETS, p)
        sort!(SUBSETS[p](db))
    elseif haskey(MATRIXCLASS, p)
        sort!(MATRIXCLASS[p])
    elseif haskey(usermatrixclass, p)
        sort!(usermatrixclass[p])
    else
        String[]
    end
end

function list(p::AbstractString, db::MatrixDatabase=MATRIX_DB)
    if startswith(p, "//")
        p = p[3:end]
        depth = 0
    else
        m = match(r"^(([^/]+/)+)(/|$)", p)
        depth = m !== nothing ? count(x == '/' for x in m.captures[1]) : -1
    end
    p = replace(p, r"//+" => '/')
    endswith(p, '/') && ( p = string(p, "**") )

    if occursin(r"[*?.]", p)
        p = replace(p, "**" => "([^/]+/)\x01[^/]\x01")
        p = replace(p, '*' => "[^/]*")
        p = replace(p, '?' => "[^/]")
        p = replace(p, '.' => "[.]")
        p = replace(p, '\x01' => '*')
    end

    r = Regex(string('^', p, '$'))
    if depth >= 0
        length(p) == 1 && ( p = ".*/" )
        listdir(r, depth, db)
    else
        list(r, db)
    end
end

list_all(db::MatrixDatabase) = sort!(collect(keys(db.data)))
function list_remote(db::MatrixDatabase)
    collect(d.name for d in values(db.data) if d isa RemoteMatrixData)
end
function list_loaded(db::MatrixDatabase)
    collect(d.name for d in values(db.data) if d isa RemoteMatrixData && !isempty(d.metadata))
end
function list_unloaded(db::MatrixDatabase)
    collect(d.name for d in values(db.data) if d isa RemoteMatrixData && isempty(d.metadata))
end
list_local(db::MatrixDatabase) = union(collect(keys(MATRIXDICT)), keys(USERMATRIXDICT))
list_builtin(db::MatrixDatabase) = collect(keys(MATRIXDICT))
list_user(db::MatrixDatabase) = collect(keys(USERMATRIXDICT))

const SUBSETS = Dict(
                     :remote => list_remote,
                     :loaded => list_loaded,
                     :unloaded => list_unloaded,
                     :local => list_local,
                     :builtin => list_builtin,
                     :user => list_user,
                     :all => list_all,
)

flatten(a) = collect(Iterators.flatten(a))
list(i::Integer, db::MatrixDatabase=MATRIX_DB) = list(Regex(string('^', aliasid(i), '$')), db)
function list(r::OrdinalRange, db::MatrixDatabase=MATRIX_DB)
    listdb(r) = list(r, db)
    sort!(flatten(listdb.(aliasid.(r) ∩ keys(db.aliases))))
end
function list(r::AbstractVector, db::MatrixDatabase=MATRIX_DB)
    listdb(r) = list(r, db)
    sort!(flatten(listdb.(r)))
end
function list(r::Tuple, db::MatrixDatabase=MATRIX_DB)
    listdb(r) = list(r, db)
    sort!(intersect(listdb.(r)...))
end

const Pattern = Union{AbstractString,Regex,OrdinalRange,Integer,Symbol,AbstractVector,Tuple}

"return information about all matrices selected by pattern"
function info(p::Pattern, db::MatrixDatabase=MATRIX_DB)
    mdbuffer = Markdown.MD([])
    md = mdbuffer
    for name in list(p)
        try
            md = info(db.data[name])
            #=
            TODO add names of meta-files to md
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
see also [`@matrix`](@ref).
"""
rhs(p::Pattern, db::MatrixDatabase=MATRIX_DB) = rhs(mdopen(p))

"""
    solution(p[, db])
return solution of problem corresponding to right hand side or `nothing`.
see also [`@matrix`](@ref).
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

###
# vintage API
###
"""
    matrixdepot(pattern, [(:search | :read | :get),])
+ search all matrix names according to `pattern` - see also [`@list`](@ref)
+ read matrix named by `pattern` - see also [`@matrix`](@ref)
+ get all matrices named according to `pattern`- see also [`@load`](@ref)
"""
function matrixdepot(p::Pattern, s::Symbol)
    if symbol in (:s, :search)
        list(p, MATRIX_DB)
    elseif s in (:g, :get)
        load(p, MATRIX_DB)
    elseif s in (:r, :read)
        matrix(p, MATRIX_DB)
    else
        throw(ArgumentError("unknown symbol $s"))
    end
end

"""
    matrixdpot(pattern)
return formatted info about all matrices with names matching `pattern`.
"""
matrixdepot(p::Pattern) = info(p, MATRIX_DB)

"""
    matrixdepot(string)
if `string` is a group name, return list of members of group

otherwise return info about pattern given by `string`
"""
function matrixdepot(p::AbstractString)
# cover the cases, a group name is given as a string and "data"
    db = MATRIX_DB
    p == "data" && (return list(:loaded, db))
    if isempty(list(p, db))
        list(Symbol(p), db)
    else
        info(p, db)
    end
end

"""
    matrixdepot(name, arg, args...)

Generate built-in or user-defined matrix named `name` with (at least one) argument    
"""
matrixdepot(p::Pattern, args...) = matrix(p, MATRIX_DB, args...)

"""
    matrixdepot()
Overview about matrices.
"""
function matrixdepot()
    #display(overview(MATRIX_DB))
    println(overview(MATRIX_DB))
end

