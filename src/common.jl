########################
# helper functions
########################
using Markdown

argerr(s::AbstractString) = throw(ArgumentError(s))
daterr(s::AbstractString) = throw(DataError(s))
parserr(s::AbstractString) = throw(Meta.ParseError(s))

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
    prop in keys(MATRIXCLASS) && daterr("$prop can not be modified.")

    user = abspath(MY_DEPOT_DIR, "group.jl")
    s = read(user, String)          # read complete file into s
    rg = Regex(repr(prop) * r"\W*=>\W*(\[.*\]\W*,\W*\n)".pattern)
    ppos = findfirst(rg, s)         # locate the prop in user.jl to remove.
    if ppos !== nothing
        start_char = first(ppos) - 1    # the start of the line
        end_char = last(ppos)           # the end of the line
    else
        ppos = findnext(r"\);", s, 1)
        start_char = ppos !== nothing ? first(ppos) - 1 : length(s)
        end_char = start_char
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

function addtogroup(dir::Dict, groupname::Symbol, f::Function)
    if groupname in keys(dir)
        fn = fname(f)
        gr = dir[groupname]
        fn in gr || push!(gr, fn)
        true
    else
        false
    end
end
function include_generator(::Type{Group}, groupname::AbstractString, f::Function)
    groupname = Symbol(groupname)
    addtogroup(MATRIXCLASS, groupname, f) ||
    addtogroup(usermatrixclass, groupname, f) ||
    argerr("$(groupname) is not a group in MatrixDepot, use
              @addgroup to add this group")
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
listdir([db,] p::AbstractString)

list directories and the number of matrices contained in them.
get an overview of the count of names down directories.
return a list with summary information for directories in matrix name space.
The input argument is split into 2 patterns by the first double slash `"//"`.
The whole string (with multiple slashes reduced to single slashes) determines
a subset of all matrix names. They are then grouped by the first pattern and
for each different group value the number of names in the subset is counted.
A final `/` is replaced by `"//**"`.

E.g.
+ `listdir("/*")`     - count names without a `/`. 
+ `listdir("/k*")`    - count names without `/` starting with `k*`.
+ `listdir("*//*")`   - count names with one directory part (uf-collection)
+ `listdir("*/*//*")` - count names with two directory parts (mm-collection)
+ `listdir("*//*/*")` - count names with two directory parts (mm-collection)
+ `listdir("Har*//*/*")` - restrict to directories starting with "Har"
+ `listdir("Har*/*//*")` - all subdirectoreis of the previous
"""
listdir(p::AbstractString) = listdir(MATRIX_DB, p)
function listdir(db::MatrixDatabase, p::AbstractString)
    r = findfirst(r"/+", p)
    if r !== nothing && first(r) == 1
        p = p[last(r)+1:end]
        depth = 0
    else
        m = match(r"^(([^/]+/)+)(/|$)", p)
        depth = m !== nothing ? count(x == '/' for x in m.captures[1]) : -1
    end
    p = replace(p, r"//+" => '/')
    endswith(p, '/') && ( p = string(p, "**") )
    r = shell_to_regex(p, false)
    if depth >= 0
        length(p) == 1 && ( p = ".*/" )
        listdir(db, r, depth)
    else
        argerr("pattern '$p' needs '//' in the middle or '/' at start or end")
    end
end
listdir(r::Regex, depth::Int) = listdir(MATRIX_DB, r)
function listdir(db::MatrixDatabase, r::Regex, depth::Int)
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

"""
    list(p::Pattern[, db]
return a vector of full matrix names where name or alias match given pattern.
`p` can be one of the following:
+ a plain string (without characters `*` and `?`) which must match exactly
+ a string containing `*` and `?` acting like a shell path pattern
+ a regular expression
+ an integer matching equivalent to the alias string `"#\$p"`
+ a range of integers
+ a group name expressed as a symbol e.g. `:local`, `:all`, `:illcond`, `posdef`
+ the name of a predicate function `f(::MatrixData)::Bool`, e.g. `symmetric`, ...
+ a vector of patterns meaning the union
+ a tuple of patterns meaning the intersection
"""
list(r::Regex) = list(MATRIX_DB, r)
function list(db::MatrixDatabase, r::Regex)
    result = AbstractString[]
    for name in keys(db.data)
        if match(r, name) !== nothing
            push!(result, name)
        end
    end
    sort!(result)
    unique!(result)
    result
end

list(p::Symbol) = list(MATRIX_DB, p)
function list(db::MatrixDatabase, p::Symbol)
    if haskey(SUBSETS, p)
        sort!(SUBSETS[p](db))
    elseif haskey(MATRIXCLASS, p)
        sort!(MATRIXCLASS[p])
    elseif haskey(usermatrixclass, p)
        sort!(usermatrixclass[p])
    else
        EMPTY_PATTERN
    end
end

"""
    shell_to_regex

return a regular expression if shell pattern characters `"*?]"` are contained in 
string, otherwise return string.
"""
function shell_to_regex(p::AbstractString, retain_pure::Bool)
    regex(p) = Regex(string('^', p, '$'))
    if occursin(r"[*?.]", p)
        p = replace(p, "**" => "\x02\x01")
        p = replace(p, '*' => "[^/]*")
        p = replace(p, '?' => "[^/]")
        p = replace(p, '.' => "[.]")
        p = replace(p, '\x01' => '*')
        p = replace(p, '\x02' => '.')
        regex(p)
    else
        retain_pure ? p : regex(p)
    end
end

list(p::AbstractString) = list(MATRIX_DB, p)
function list(db::MatrixDatabase, p::AbstractString)
    p = replace(p, r"//+" => '/')
    r = shell_to_regex(p, true)
    r isa Regex ? list(db, r) : haskey(db.data, r) ? [r] : EMPTY_PATTERN
end


list(p::Alias) = list(MATRIX_DB, p)
list(db::MatrixDatabase, p::Alias) = list(aliasresolve(db, p))

list(p::Not) = list(MATRIX_DB, p)
list(db::MatrixDatabase, p::Not) = setdiff(list(db,()), list(db, p.pattern)) 

list(r::AbstractVector) = list(MATRIX_DB, r)
function list(db::MatrixDatabase, r::AbstractVector)
    listdb(r) = list(db, r)
    unique!(sort!(flatten(listdb.(r))))
end
list(::Tuple{}) = list_all(MATRIX_DB)
list(db::MatrixDatabase, ::Tuple{}) = list_all(db)
list(r::Tuple) = list(MATRIX_DB, r)
function list(db::MatrixDatabase, r::Tuple)
    y, st = iterate(r)
    res = list(y)
    while !isempty(res) && (x = iterate(r, st)) !== nothing
        y, st = x
        intersect!(res, list(y))
    end
    res
end
list(pred::Function) = list(MATRIX_DB, pred)
function list(db::MatrixDatabase, pred::Function)
    dali = [ data for data in values(db.data) if pred(data) ]
    sort!(getfield.(dali, :name))
end

## internal list special cases
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

"return information about all matrices selected by pattern"
info(p::Pattern) = info(MATRIX_DB, p)
function info(db::MatrixDatabase, p::Pattern)
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

function verify_loaded(db::MatrixDatabase, data::RemoteMatrixData)
    if isempty(data.metadata)
        loadmatrix(db, data)
    end
    data
end
verify_loaded(db::MatrixDatabase, data::MatrixData) = data
mdatav(db::MatrixDatabase, p::Pattern) = verify_loaded(db, mdata(db, p))

"""
    matrix(p, [db,] args...)
return a matrix given a name
* `p` is a unique problem identifier (`String`, `Regex`, `Integer`)
* `db` is optional `MatrixDatabase` defaulting to global `MATRIX_DB`
* `args` are optional arguments - only used for generated matrices
"""
matrix(p::Pattern, args...) = matrix(MATRIX_DB, p, args...)
matrix(db::MatrixDatabase, p::Pattern, args...) = matrix(mdatav(db, p), args...)

"""
    rhs(p[, db]
return right hand side of problem or `nothing`.
see also [`@matrix`](@ref).
"""
rhs(p::Pattern) = rhs(MATRIX_DB, p)
rhs(db::MatrixDatabase, p::Pattern) = rhs(mdatav(db, p))

"""
    solution(p[, db])
return solution of problem corresponding to right hand side or `nothing`.
see also [`@matrix`](@ref).
"""
solution(p::Pattern) = solution(MATRIX_DB, p)
solution(db::MatrixDatabase, p::Pattern) = solution(mdatav(db, p))

"""
    metareader(p, metafile::String[, db])
read metadata identified by name `metafile`.
"""
metareader(p::Pattern, metafile::AbstractString) = metareader(MATRIX_DB, p, metafile)
function metareader(db::MatrixDatabase, p::Pattern, metafile::AbstractString)
    metareader(mdatav(db, p), metafile)
end

"""
    load(pattern[, db])
load data from remote repository for all problems matching pattern"
return the number of successfully loaded matrices.
"""
load(p::Pattern) = load(MATRIX_DB, p)
function load(db::MatrixDatabase, p::Pattern)
    n = 0
    for name in list(p)
        try
            loadmatrix(db, db.data[name])
            n += 1
        catch ex
            @warn "could not load $name: $ex"
        end
    end
    n
end

"""
    mdata(pattern[ ,db])

make sure that data files are loaded-
return `MatrixData` object, which can be used with data access functions.
The data cache is activated for `RemoteMatrixData`. see [`@mdclose`](@ref).
If the pattern has not a unique resolution, an error is thrown.
"""
mdopen(p::Pattern; cache::Bool=false) = mdopen(MATRIX_DB, p, cache=cache)
function mdopen(db::MatrixDatabase, p::Pattern; cache::Bool=false)
    mdopen(mdatav(db, p), cache=cache)
end

"""
    mdata(db, pattern)
return unique `MatrixData` object according to pattern.
"""
mdata(p::Pattern) = mdata(MATRIX_DB, p)
function mdata(db::MatrixDatabase, p::Pattern)
    li = list(p)
    length(li) == 0 && daterr("no matrix according to $p found")
    length(li) > 1  && daterr("pattern not unique: $p -> $li")
    db.data[li[1]]
end

"""
    metadata(p::Pattern[, db])
return copy of metadata if pattern is unique
"""
metadata(data::RemoteMatrixData) = copy(data.metadata)
metadata(data::MatrixData) = String[]
metadata(p::Pattern) = metadata(MATRIX_DB, p)
metadata(db::MatrixDatabase, p::Pattern) = metadata(mdata(db, p))

"""
    mdopen(data::MatrixData)
Enable data caching for `RemoteMatrixData` and return data.
"""
mdopen(data::RemoteMatrixData; cache::Bool=false) = (data.status[] = cache; data)
mdopen(data::MatrixData; cache::Bool=false) = data

"""
    mdclose(data::MatrixData)
clean cached data and disable caching - return data.
"""
function mdclose(data::RemoteMatrixData)
    data.status[] = false
    empty!(data.datacache)
    data
end
mdclose(data::MatrixData) = data
mdclose(p::Pattern) = mdclose(MATRIX_DB, p)
function mdclose(db::MatrixDatabase, p::Pattern)
    mdatadb(name) = mdata(db, name)
    mdclose.(mdatadb.(list(db, p)))
end

###
# vintage API
###

# convert input string to shell pattern in certain cases
function convert_pattern(p::AbstractString)
    occursin(r"[]*?/]", p) ? p : ["**/$p", p]
end

"""
    matrixdepot(pattern, [(:search | :read | :get),])
+ search all matrix names according to `pattern` - see also [`@list`](@ref)
+ read matrix named by `pattern` - see also [`@matrix`](@ref)
+ get all matrices named according to `pattern`- see also [`@load`](@ref)
"""
function matrixdepot(p::AbstractString, s::Symbol; meta::Bool=false)
    p1 = convert_pattern(p)
    if s in (:s, :search)
        list(MATRIX_DB, p1)
    elseif s in (:g, :get)
        load(MATRIX_DB, p1)
    elseif s in (:r, :read)
        data = mdata(MATRIX_DB, p)
        !meta ? matrix(data) : meta_names(data)
    else
        argerr("unknown symbol $s")
    end
end

"""
    meta_names(data::MatrixData)
return dictionary of file names associated with metadata names.
"""
function meta_names(data::RemoteMatrixData)
    mlist = data.metadata
    path(name) = joinpath(dirname(matrixfile(data)), name)
    abbr(name) = rsplit(name, '.'; limit=2)[1]
    Dict{String,String}([ abbr(x) => path(x) for x in mlist])
end
meta_names(data::MatrixData) = Dict{String,String}()

"""
    matrixdepot(string)
if `string` is a group name, return list of members of group
The string "all" returns all generated (=local) matrices.
The string "data" returns all loaded remote matrices.

otherwise return info about pattern given by `string`.
"""
function matrixdepot(p::AbstractString)
# cover the cases, a group name is given as a string and "data"
    db = MATRIX_DB
    p == "data" && (return list(db, :loaded))
    p == "all" && (return list(db, :local))
    p1 = convert_pattern(p)
    if isempty(list(db, p1))
        p1= replace(p, '-' => "")
        res = list(db, Symbol(p1))
        isempty(res) ? daterr("not found '$p'") : res
    else
        info(db, p1)
    end
end

"""
    matrixdepot(n:Integer)
return name of generated matrix numbered by n.
hint: the list API returns an array about the remote matric with id n.
"""
function matrixdepot(n::Integer)
    locli = list(:local)
    1 <= n <= length(locli) || daterr("1 <= $n <= $(length(locli)) required.")
    locli[n]
end
"""
    matrixdepot(range or tuple of integer or range arguments)
return array of names corresponding to the local matrices.
"""
function matrixdepot(r::AbstractRange{<:Integer})
    matrixdepot.(r)
end
function matrixdepot(n::Union{Integer,AbstractRange{<:Integer}},
                     n2::Union{Integer,AbstractRange{<:Integer}}...)

    [matrixdepot(n); matrixdepot(n2...)]
end

"""
    matrixdepot(name, arg, args...)

Generate built-in or user-defined matrix named `name` with (at least one) argument    
"""
matrixdepot(p::Pattern, args...) = matrix(MATRIX_DB, p, args...)

"""
    matrixdepot(group::String...)
return intersection of output of individual group names.
"""
function matrixdepot(g::AbstractString, h::AbstractString...)
    first = matrixdepot(g)
    first = first isa Vector{<:AbstractString} ? first : [g]
    second = matrixdepot(h...)
    second = second isa Vector{<:AbstractString} ? second : [h...]
    intersect(first, second)
end

"""
    matrixdepot()
Overview about matrices.
"""
function matrixdepot()
    #display(overview(MATRIX_DB))
    println(overview(MATRIX_DB))
end
