########################
# helper functions
########################

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
    append!(groups, sort!(collect(keys(SUBSETS))))
    append!(groups, sort!(collect(keys(MATRIXCLASS))))
    append!(groups, sort!(collect(keys(usermatrixclass))))
    groups
end

#######################
# matrix group
#######################

# write one property association
function propline(io::IO, propname, matnames)
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
function include_generator(::Type{Group}, groupname::Symbol, f::Function)
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
listdir(p::AbstractString, xp::Pattern=()) = listdir(MATRIX_DB, p, xp)
function listdir(db::MatrixDatabase, p::AbstractString, xp::Pattern)
    check_symbols(xp)
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
        listdir(db, r & xp, depth)
    else
        argerr("pattern '$p' needs '//' in the middle or '/' at start or end")
    end
end
function listdir(db::MatrixDatabase, r::Pattern, depth::Int)
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
    listdata([db,] p::Pattern)
Return an array of `MatrixData` objects according to matched patterns.
"""
listdata(p::Pattern) = listdata(MATRIX_DB, p)
listdata(db::MatrixDatabase, p::Pattern) = mdata.(MatrixDepot.list(db, p))

"""
    list([db,] p::Pattern)
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
list(p::Pattern) = list(MATRIX_DB, p)
is_all(res::Vector) = length(res) == 1 && res[1] == ""

function list(db::MatrixDatabase, p::Pattern)
    res = list!(db, [""], p)
    is_all(res) ? list_all(db) : res
end
function list!(db::MatrixDatabase, res::Vector{String}, r::Regex)
    isempty(res) && return res
    if is_all(res)
        empty!(res)
        for name in keys(db.data)
            if match(r, name) !== nothing
                push!(res, name)
            end
        end
        sort!(res)
    else
        for i = 1:length(res)
            name = res[i]
            if match(r, name) === nothing
                res[i] = ""
            end
        end
        for i = length(res):-1:1
            if res[i] == ""
                deleteat!(res, i)
            end
        end
    end
    res
end

function list!(db::MatrixDatabase, res::Vector{String}, p::Symbol)
    isempty(res) && return res
    x = if haskey(SUBSETS, p)
        SUBSETS[p](db)
    elseif haskey(MATRIXCLASS, p)
        MATRIXCLASS[p]
    elseif haskey(usermatrixclass, p)
        usermatrixclass[p]
    else
        argerr("unknown group name '$p'")
        # EMPTY_PATTERN
    end
    if is_all(res)
        empty!(res)
        append!(res, sort(x))
    else
        intersect!(res, x)
    end
end

"""
    shell_to_regex

return a regular expression if shell pattern characters `"*?]"` are contained in
string, otherwise return string.
If no `'/'` is contained in p and p is not "*", insert `"(.*/)?"` in regular expression.
"""
function shell_to_regex(p::AbstractString, retain_pure::Bool)
    regex(p) = Regex(string('^', p, '$'))
    p = replace(p, r"//+" => '/')
    # p = p == "*" || '/' in p ? p : string("(**/)\x03", p)
    if occursin(r"[*?.]", p)
        p = replace(p, "**" => "\x01\x02")
        p = replace(p, '*' => "[^/]*")
        p = replace(p, '?' => "[^/]")
        p = replace(p, '.' => "[.]")
        p = replace(p, '\x01' => '.')
        p = replace(p, '\x02' => '*')
        p = replace(p, '\x03' => '?')
        regex(p)
    else
        retain_pure ? p : regex(p)
    end
end

function singlist!(db::MatrixDatabase, res::Vector{String}, p::AbstractString)
    if is_all(res)
        if haskey(db.data, p)
            res[1] = p
            res
        else
            empty!(res)
        end
    else
        if p in res
            empty!(res)
            push!(res, p)
        else
            empty!(res)
        end
    end
end

function list!(db::MatrixDatabase, res::Vector{String}, p::AbstractString)
    isempty(res) && return res
    r = shell_to_regex(p, true)
    r isa Regex ? list!(db, res, r) : singlist!(db, res, r)
end

list!(db::MatrixDatabase, res::Vector{String}, p::Alias) = list!(db, res, aliasresolve(db, p))

# If res is symbolically [""], populate it with the set of all available names
function resall!(db::MatrixDatabase, res::Vector{String})
    if is_all(res)
        x = list_all(db)
        resize!(res, length(x))
        copyto!(res, x)
    end
    res
end

function list!(db::MatrixDatabase, res::Vector{String}, p::Not)
    isempty(res) && return res
    cres = copy(res)
    resall!(db, res)
    setdiff!(res, list!(db, cres, p.pattern))
end

# logical OR
function list!(db::MatrixDatabase, res::Vector{String}, r::AbstractVector)
    isempty(r) && return empty!(res)
    check_symbols(r)
    isempty(res) && return res
    length(r) == 1 && return list!(db, res, r[1])
    cres = copy(res)
    list!(db, res, r[1])
    for y in r[2:end]
        union!(res, list!(db, copy(cres), y))
    end
    sort!(res)
end

list!(db::MatrixDatabase, res::Vector{String}, ::Tuple{}) = res

# logical AND 
function list!(db::MatrixDatabase, res::Vector{String}, r::Tuple)
    isempty(res) && return res
    check_symbols(r)
    for y in sort_by_type(r)
        list!(db, res, y)
        isempty(res) && break
    end
    res
end

function list!(db::MatrixDatabase, res::Vector{String}, pred::Function)
    isempty(res) && return res
    resall!(db, res)
    filter!(k -> pred(db.data[k]), res)
end

# return a vector with re-arranged the contents of itr
# the elements of its must be of type Pattern
function sort_by_type(itr)
    f(p::AbstractString) = true
    f(p::Symbol) = true
    f(p::AbstractVector) = all(f.(p))
    f(p::Tuple) = all(f.(p))
    f(p::Pattern) = false
    vall = collect(itr)
    vend = filter(!f, vall)
    filter!(f, vall)
    append!(vall, vend)
    vall
end

## internal list special cases
list_all(db::MatrixDatabase) = sort!(String.(collect(keys(db.data))))
list_local(db::MatrixDatabase) = union(collect(keys(MATRIXDICT)), keys(USERMATRIXDICT))
list_builtin(db::MatrixDatabase) = collect(keys(MATRIXDICT))
list_user(db::MatrixDatabase) = collect(keys(USERMATRIXDICT))

const SUBSETS = Dict(
                     :local => list_local,
                     :builtin => list_builtin,
                     :user => list_user,
                     :all => list_all,
)

function verify_loaded(db::MatrixDatabase, data::RemoteMatrixData)
    if isempty(data.metadata) || !isfile(matrixfile(data))
        loadmatrix(data)
    end
    data
end
verify_loaded(db::MatrixDatabase, data::MatrixData) = data

function verify_loadinfo(data::RemoteMatrixData)
    file = matrixfile(data)
    if !isfile(file)
        file = matrixinfofile(data)
        if !isfile(file)
            loadinfo(data)
        end
    end
    file
end

mdatav(db::MatrixDatabase, p::Pattern) = verify_loaded(db, mdata(db, p))

"""
    load(pattern[, db])

Load data from remote repository for all problems matching pattern.

Return the number of successfully loaded matrices. 
"""
load(p::Pattern) = load(MATRIX_DB, p)
function load(db::MatrixDatabase, p::Pattern)
    check_symbols(p)
    n = 0
    for name in list(p)
        try
            n += loadmatrix(db.data[name])
        catch ex
            ex isa InterruptException && rethrow()
            @warn "could not load $name: $ex"
        end
    end
    n
end

"""
    mdopen([db,] pattern)
    mdopen(f, [db,] pattern)

Return `MatrixDescriptor` object, which can be used with data access functions.

Make sure that data files are loaded.
Keeps a cache of already delivered matrices and metadata.
If the pattern has not a unique resolution, an error is thrown.
"""
mdopen(p::Pattern, args...) = mdopen(MATRIX_DB, p, args...)
function mdopen(db::MatrixDatabase, p::Pattern, args...)
    _mdopen(mdatav(db, p), args...)
end

mdopen(f::Function, p::Pattern, args...) = mdopen(f, MATRIX_DB, p, args...)
function mdopen(f::Function, db::MatrixDatabase, p::Pattern, args...)
    data = _mdopen(mdatav(db, p), args...)
    f(data)
end

"""
    mdata(db, pattern)

Return unique `MatrixData` object according to pattern.
"""
mdata(p::Pattern) = mdata(MATRIX_DB, p)
function mdata(db::MatrixDatabase, p::Pattern)
    check_symbols(p)
    li = list(db, p)
    length(li) == 0 && daterr("no matrix according to $p found")
    length(li) > 1  && daterr("pattern not unique: $p -> $li")
    db.data[li[1]]
end

"""
    metadata([db, ], p::Pattern)

Return copy of list of metadata names. Pattern must be unique.
"""
metadata(mdesc::MatrixDescriptor) = metadata(mdesc.data)
metadata(data::RemoteMatrixData) = copy(data.metadata)
metadata(data::MatrixData) = String[]
metadata(p::Pattern) = metadata(MATRIX_DB, p)
metadata(db::MatrixDatabase, p::Pattern) = metadata(mdata(db, p))

_mdopen(data::RemoteMatrixData)= MatrixDescriptor(data)
function _mdopen(data::GeneratedMatrixData, args...)
    # verify_callable(data.func, args...)
    md = MatrixDescriptor(data, args...)
    md.A
    md
end

verify_callable(f::Function, args...) = verify_callable(f, Float64, args...)
function verify_callable(f::Function, t::Type, args...)
    ml = methods(f, typeof.((t, args...)))
    if isempty(ml) || verify_types(ml.ms[1].sig)
        throw(MethodError(f, (t, args...)))
    end
end

verify_types(sig::UnionAll) = verify_types(sig.body)
function verify_types(sig::DataType)
    types = sig.types
    length(types) == 3 && types[3] == Vararg{Any} && types[2] isa Type
end

###
# convenience API
###

"""
    matrixdepot(p::Pattern, args...)

Return matrix according to pattern or local matrix according to name and arguments.

If not loaded, load remote matrix first.
`p` must be a unique pattern (match only one name). The presence of arguments makes
sense only if the pattern matches the name of a generated (=local) matrix.

Only the matrix part is delivered, also in the local cases, where the underlying
function returns a structure containing matrix and vectors.
Use `md = mdopen; md.A, md.b ...`
to access those objects.
"""
matrixdepot(p::Pattern, args...) = matrixdepot(MATRIX_DB, p, args...)
function matrixdepot(db::MatrixDatabase, p::Pattern, args...)
    mdopen(db, p, args...) do md
        md.A
    end
end
