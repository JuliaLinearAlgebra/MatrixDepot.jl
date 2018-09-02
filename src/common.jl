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

function include_generator(::Type{Group}, groupname::AbstractString, f::Function)
    groupname = Symbol(groupname)
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
+ a group name expressed as a symbol e.g. `:local`, `:all`, `:illcond`, `posdef`
+ the name of a predicate function `f(::MatrixData)::Bool`, e.g. `symmetric`, ...
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
list(i::Integer, db::MatrixDatabase=MATRIX_DB) = list(Regex(string('^', aliasname(i), '$')), db)
function list(r::OrdinalRange, db::MatrixDatabase=MATRIX_DB)
    listdb(r) = list(r, db)
    unique!(sort!(flatten(listdb.(aliasname.(r) ∩ keys(db.aliases)))))
end
function list(r::AbstractVector, db::MatrixDatabase=MATRIX_DB)
    listdb(r) = list(r, db)
    unique!(sort!(flatten(listdb.(r))))
end
function list(r::Tuple, db::MatrixDatabase=MATRIX_DB)
    listdb(r) = list(r, db)
    sort!(intersect(listdb.(r)...))
end
function list(pred::Function, db::MatrixDatabase=MATRIX_DB)
    dali = [ data for data in values(db.data) if pred(data) ]
    sort!(getfield.(dali, :name))
end

const Pattern = Union{AbstractString,Regex,OrdinalRange,Integer,Symbol,
                      AbstractVector,Tuple,Function}

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

function verify_loaded(data::RemoteMatrixData, db::MatrixDatabase)
    if isempty(data.metadata)
        loadmatrix(data, db)
    end
    data
end
verify_loaded(data::MatrixData, db::MatrixDatabase) = data
mdatav(p::Pattern, db::MatrixDatabase) = verify_loaded(mdata(p, db), db)

"""
    matrix(p, [db,] args...)
return a matrix given a name
* `p` is a unique problem identifier (`String`, `Regex`, `Integer`)
* `db` is optional `MatrixDatabase` defaulting to global `MATRIX_DB`
* `args` are optional arguments - only used for generated matrices
"""
matrix(p::Pattern, db::MatrixDatabase=MATRIX_DB, args...) = matrix(mdatav(p, db), args...)
matrix(p::Pattern, args...) = matrix(mdatav(p, MATRIX_DB), args...)

"""
    rhs(p[, db]
return right hand side of problem or `nothing`.
see also [`@matrix`](@ref).
"""
rhs(p::Pattern, db::MatrixDatabase=MATRIX_DB) = rhs(mdatav(p, db))

"""
    solution(p[, db])
return solution of problem corresponding to right hand side or `nothing`.
see also [`@matrix`](@ref).
"""
solution(p::Pattern, db::MatrixDatabase=MATRIX_DB) = solution(mdatav(p, db))

"""
    metareader(p, metafile::String[, db])
read metadata identified by name `metafile`.
"""
function metareader(p::Pattern, metafile::AbstractString, db::MatrixDatabase=MATRIX_DB)
    metareader(mdatav(p, db), metafile)
end

"""
    load(pattern[, db])
load data from remote repository for all problems matching pattern"
return the number of successfully loaded matrices.
"""
function load(p::Pattern, db::MatrixDatabase=MATRIX_DB)
    n = 0
    for name in list(p)
        try
            loadmatrix(db.data[name], db)
            n += 1
        catch ex
            @warn "could not load $name: $ex"
        end
    end
    n
end

"""
    mdata(pattern[ ,db])
return `MatrixData` object, which can be used with data access functions.
The data cache is activated for `RemoteMatrixData`. see [`@mdclose`](@ref).
If the pattern has not a unique resolution, an error is thrown.
"""
function mdopen(p::Pattern, db::MatrixDatabase=MATRIX_DB)
    mdopen(mdatav(p, db)) # set status flag
end

"""
    mdata(pattern, db)
return unique `MatrixData` object according to pattern.
"""
function mdata(p::Pattern, db::MatrixDatabase)
    li = list(p)
    length(li) == 0 && daterr("no matrix according to $p found")
    length(li) > 1  && daterr("pattern not unique: $p -> $li")
    db.data[li[1]]
end

"""
    metadata(p::Pattern[, db])
return copy of metadata if pattern is unique
"""
metadata(p::Pattern, db::MatrixDatabase=MATRIX_DB) = copy(mdata(p, db).metadata)

"""
    mdopen(data::MatrixData)
Enable data caching for `RemoteMatrixData` and return data.
"""
mdopen(data::RemoteMatrixData) = (data.status[] = true; data)
mdopen(data::MatrixData) = data

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

###
# vintage API
###

# convert input string to shell pattern in certain cases
function convert_pattern(p::AbstractString)
    list(occursin(r"[]*?/]", p) ? p : ["**/$p", p], MATRIX_DB)
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
        list(p1, MATRIX_DB)
    elseif s in (:g, :get)
        load(p1, MATRIX_DB)
    elseif s in (:r, :read)
        data = mdata(p, MATRIX_DB)
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
    p == "data" && (return list(:loaded, db))
    p == "all" && (return list(:local, db))
    p1 = convert_pattern(p)
    if isempty(list(p1, db))
        p1= replace(p, '-' => "")
        res = list(Symbol(p1), db)
        isempty(res) ? daterr("not found '$p'") : res
    else
        info(p1, db)
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
matrixdepot(p::Pattern, args...) = matrix(p, MATRIX_DB, args...)

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

###
# Syntactic sugar
###
import Base: &, |, !
(&)(p::Pattern, q::Pattern...) = tuple(p, q...)
(|)(p::Pattern, q::Pattern...) = vcat(p, q...)

###
# Predefined predicates for MatrixData
###
export isgeneral, issymmetric, isskew, ishermitian
export iscomplex, isreal, isinteger, ispattern
export isremote, islocal, isloaded, isunloaded, isbuiltin, isuser
export predm, predn, prednz, predmn

import Base: isreal, isinteger
import LinearAlgebra: issymmetric, ishermitian

function issymmetry(data::RemoteMatrixData, T::Type{<:MMSymmetry})
    data.properties[] !== nothing && data.properties[].symmetry isa T
end
function isfield(data::RemoteMatrixData, T::Type{<:MMField})
    data.properties[] !== nothing && data.properties[].field isa T
end

isgeneral(data::RemoteMatrixData) = issymmetry(data, MMSymmetryGeneral)
issymmetric(data::RemoteMatrixData) = issymmetry(data, MMSymmetrySymmetric)
isskew(data::RemoteMatrixData) = issymmetry(data, MMSymmetrySkewSymmetric)
ishermitian(data::RemoteMatrixData) = issymmetry(data, MMSymmetryHermitian)
isgeneral(data::MatrixData) = !issymmetric(data) && !isskew(data) && !ishermitian(data)
issymmetric(data::MatrixData) = data.name in list(:symmetric)
isskew(data::MatrixData) = false
ishermitian(data::MatrixData) = false

iscomplex(data::RemoteMatrixData) = isfield(data, MMFieldComplex)
isreal(data::RemoteMatrixData) = isfield(data, MMFieldReal)
isinteger(data::RemoteMatrixData) = isfield(data, MMFieldInteger)
ispattern(data::RemoteMatrixData) = isfield(data, MMFieldPattern)
iscomplex(data::MatrixData) = false
isreal(data::MatrixData) = false
isinteger(data::MatrixData) = false
ispattern(data::MatrixData) = false

hasproperties(data::RemoteMatrixData) = data.properties[] !== nothing
hasproperties(data::MatrixData) = false
isremote(data::RemoteMatrixData) = true
isremote(data::MatrixData) = false
isloaded(data::RemoteMatrixData) = !isempty(data.metadata)
isloaded(data::MatrixData) = false
isunloaded(data::RemoteMatrixData) = isempty(data.metadata)
isunloaded(data::MatrixData) = false
isuser(data::GeneratedUserMatrixData) = true
isuser(data::MatrixData) = false
isbuiltin(data::GeneratedBuiltinMatrixData) = true
isbuiltin(data::MatrixData) = false
islocal(data::GeneratedMatrixData) = true
islocal(data::MatrixData) = false

predm(f::Function) = data::MatrixData -> hasproperties(data) && f(data.properties[].m)
predn(f::Function) = data::MatrixData -> hasproperties(data) && f(data.properties[].n)
prednz(f::Function) = data::MatrixData -> hasproperties(data) && f(data.properties[].nz)
function predmn(f::Function)
    data::MatrixData -> hasproperties(data) && f(data.properties[].m, data.properties[].n)
end
