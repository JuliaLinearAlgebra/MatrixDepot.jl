
# functions related to the logical operations on patterns
# and predicate functions.
#
export builtin, user, uf, tamu, mm, ¬, logical

export isgeneral, issymmetric, isskew, ishermitian
export iscomplex, isreal, isinteger, ispattern
export isremote, islocal, isloaded, isunloaded, isbuiltin, isuser
export @pred

import Base: isreal, isinteger
import LinearAlgebra: issymmetric, ishermitian

import Base: &, |, *

"""
    Pattern syntactic sugar

the logical operators `|`, `&`, and `¬` (`\\neg`)
 may be applied to all kind of `Pattern`s with the usual meaning.
+ `¬` : unary negation operator - pattern does not match (highest priority)
+ `&` : binary logical and - both patterns match
+ `|` : binary logical or - any of the pattern match (lowest priority)

+ parentheses can be used to overrule operator precedence.

+ `[p...]` is the same as `p[1] | p[2] ...`

+ `(p...)` is the same as `p[1] & p[2] ...`

+ `¬(p...) === ¬((p...))` - that is `¬(p[1] & p[2] ...)`

+ Precedence of '*' is higher that `¬` for character and string objects:
 so `¬ "a" * "b"  === ¬("a" * "b") === ¬"ab"`  also `¬'a'^2*"b" === ¬"aab"`
"""
function logical end

(&)(p::Tuple, q::Tuple) = tuple(p..., q...)
(&)(p::Tuple, q::Pattern...) = tuple(p..., q...)
(&)(p::Pattern, q::Tuple) = tuple(p, q...)
(&)(p::Pattern, q::Pattern...) = tuple(p, q...)
(|)(p::Pattern, q::Pattern...) = vcat(p, q...)
(¬)(p::Pattern) = Not(p)
(¬)(c::AbstractChar) = Not(string(c))
(¬)(p::Not) = p.pattern
Not(p::Not) = p.pattern
(*)(a::Not{<:AbstractString}, b::Union{AbstractString,AbstractChar}) = Not(a.pattern * b)
(¬)(p::Pattern...) = Not(tuple(p...))
(¬)() = EMPTY_PATTERN
(¬)(p::Vector) = length(p) == 0 ? ALL_PATTERN : Not(p)

const EMPTY_PATTERN = []
const ALL_PATTERN = ()

###
# simplification of patterns
###
is_pure_string(p::AbstractString) = !occursin(r"[]*?]", p)
is_pure_string(p::Pattern) = false
is_pure_vector(p::Vector) = all(is_pure_string.(p))
is_pure_vector(p::Pattern) = false

function simplify(v::Vector)
    w = simplify.(v)
    x = Pattern[]
    for wi in w
        wi == () && return wi
        wi isa Vector ? append!(x, wi) : push!(x, wi)
    end
    isempty(x) && return x
    unique!(x)
    pstr = filter(is_pure_string, x)
    if length(pstr) == length(x)
        length(pstr) == 1 ? pstr[1] : pstr
    elseif isempty(pstr)
        x
    else
        filter!(!is_pure_string, x)
        insert!(x, 1, prednames(pstr))
    end
end

function simplify(v::Tuple{Vararg{<:Pattern}})
    w = simplify.(v)
    x = Pattern[]
    for wi in w
        wi == Pattern[] && return wi
        wi isa Tuple ? append!(x, wi) : push!(x, wi)
    end
    isempty(x) && return ()
    unique!(x)
    pstr = filter(is_pure_string, x)
    if length(pstr) == length(x)
        length(pstr) == 1 ? pstr[1] : Pattern[]
    elseif isempty(pstr)
        x = simplify2(x)
        length(x) == 1 ? x[1] : tuple(x...)
    else
        if length(pstr) == 1
            filter!(!is_pure_string, x)
            y = filter(is_pure_vector, x)
            if all(in.(pstr[1], y))
                filter!(!is_pure_vector, x)
                insert!(x, 1, pstr[1])
                x = simplify2(x)
                length(x) == 1 ? x[1] : tuple(x...)
            else
                Pattern[]
            end
        else
            Pattern[]
        end
    end
end

function simplify2(x::Vector{<:Pattern})
    y = filter(is_pure_vector, x)
    filter!(!is_pure_vector, x)
    y = intersect(y...)
    n = length(y)
    n == 0 ? Pattern[] : n == 1 ? insert!(x, 1, y[1]) : insert!(x, 1, y)
end

function Not(pp::Function)
    if (:f,) == propertynames(pp)
        pp.f
    else
        !pp
    end
end

function simplify(p::Function)
    if (:f,) == propertynames(p) && (:f,) == propertynames(p.f)
        simplify(p.f.f)
    else
        p
    end
end

simplify(p::Pattern) = p

#=
function list(p::Tuple)
    p = simplify(p)
    k = findfirst(is_pure_vector, p)
    all = k === nothing ? keys(db.data) : p[k]
end
=#

###
# Alias resolution
###

function aliasresolve(db::MatrixDatabase, k::AbstractString)
    haskey(db.aliases, k) ? [db.aliases[k]] : String[]
end
function aliasresolve(db::MatrixDatabase, a::Alias{T,<:Integer}) where T
    aliasresolve(db, aliasname(a))
end
function aliasresolve(db::MatrixDatabase, a::Alias{T,<:AbstractVector{<:IntOrVec}}) where T
    aliasr2(x) = aliasresolve(db, x)
    collect(Iterators.flatten(aliasr2.(aliasname(a))))
end
aliasresolve(db::MatrixDatabase, a::Alias{RemoteMatrixData{TURemoteType},Colon}, ) = "*/*"
aliasresolve(db::MatrixDatabase, a::Alias{RemoteMatrixData{MMRemoteType},Colon}) = "*/*/*"
aliasresolve(db::MatrixDatabase, a::Alias{GeneratedMatrixData{:B},Colon}) = :builtin
aliasresolve(db::MatrixDatabase, a::Alias{GeneratedMatrixData{:U},Colon}) = :user

###
# Predefined predicates for MatrixData
###

builtin(p...) = Alias{GeneratedMatrixData{:B}}(p...)
user(p...) = Alias{GeneratedMatrixData{:U}}(p...)
uf(p...) = Alias{RemoteMatrixData{TURemoteType}}(p...)
tamu(p...) = Alias{RemoteMatrixData{TURemoteType}}(p...)
mm(p...) = Alias{RemoteMatrixData{MMRemoteType}}(p...)

function _issymmetry(data::RemoteMatrixData, T::Type{<:MMSymmetry})
    data.properties[] !== nothing && data.properties[].symmetry isa T
end
function _isfield(data::RemoteMatrixData, T::Type{<:MMField})
    data.properties[] !== nothing && data.properties[].field isa T
end

isgeneral(data::RemoteMatrixData) = _issymmetry(data, MMSymmetryGeneral)
issymmetric(data::RemoteMatrixData) = _issymmetry(data, MMSymmetrySymmetric)
isskew(data::RemoteMatrixData) = _issymmetry(data, MMSymmetrySkewSymmetric)
ishermitian(data::RemoteMatrixData) = _issymmetry(data, MMSymmetryHermitian)
isgeneral(data::MatrixData) = !issymmetric(data) && !isskew(data) && !ishermitian(data)
issymmetric(data::MatrixData) = data.name in list(:symmetric)
isskew(data::MatrixData) = false
ishermitian(data::MatrixData) = false

iscomplex(data::RemoteMatrixData) = _isfield(data, MMFieldComplex)
isreal(data::RemoteMatrixData) = _isfield(data, MMFieldReal)
isinteger(data::RemoteMatrixData) = _isfield(data, MMFieldInteger)
ispattern(data::RemoteMatrixData) = _isfield(data, MMFieldPattern)
iscomplex(data::MatrixData) = false
isreal(data::MatrixData) = false
isinteger(data::MatrixData) = false
ispattern(data::MatrixData) = false

hasinfo(data::RemoteMatrixData) = data.header.m > 0 && data.header.n > 0
hasinfo(data::MatrixData) = false
isremote(data::RemoteMatrixData) = true
isremote(data::MatrixData) = false
isloaded(data::RemoteMatrixData) = !isempty(data.metadata)
isloaded(data::MatrixData) = false
isunloaded(data::RemoteMatrixData) = isempty(data.metadata)
isunloaded(data::MatrixData) = false
isuser(data::GeneratedMatrixData{:U}) = true
isuser(data::MatrixData) = false
isbuiltin(data::GeneratedMatrixData{:B}) = true
isbuiltin(data::MatrixData) = false
islocal(data::GeneratedMatrixData) = true
islocal(data::MatrixData) = false

function pred(f::Function, s::Symbol...)
    data::MatrixData -> hasinfo(data) &&
    s ⊆ propertynames(data) &&
    f(getproperty.(Ref(data), s)...)
end

function prednzdev(dev::AbstractFloat=0.1)
    function f(data::RemoteMatrixData)
        n1, n2 = extremnnz(data)
        n1 -= Int(floor(n1 * dev))
        n2 += Int(floor(n2 * dev))
        isloaded(data) && ! ( n1 <= nz_num(data) <= n2 )
    end
    f(::MatrixData) = false
    f
end

"""
    check_symbols(p::Pattern)

throw `ArgumentError` if pattern uses unknown symbol as a group name.
"""
function check_symbols(p::Pattern)
    s = setdiff(filter(x->x isa Symbol, flatten_pattern(p)), group_list())
    isempty(s) || argerr("The following symbols are no group names: $s")
end

# Predicate generating macros

# extract all symbols from an expression
function extract_symbols(ex)
    s = Set{Symbol}()
    append_symbols!(::Any) = s
    append_symbols!(ex::Symbol) = push!(s, ex)
    function append_symbols!(ex::Expr)
        append_symbols!(ex.head)
        append_symbols!.(ex.args)
        s
    end
    collect(append_symbols!(ex))
end

# construct a function definition from a list of symbols and expression
function make_func(sli::AbstractVector{Symbol}, ex)
    res = :( () -> $ex )
    append!(res.args[1].args, sli)
    res
end

function make_pred(ex)
    syms = extract_symbols(ex) ∩ (:m, :n, :nnz, :dnz, :name, :id, :title, :author, :ed, :fields, :notes, :date, :kind, :metadata)
    :(pred($(make_func(syms, ex)), $(QuoteNode.(syms)...)))
end

"""
    @pred(expression)

Generate a predicate function using the expression as function body. Variable names
within the expression, which are properties of `RemoteMatrixData` (e.g. `title`, `m`, `nnz`)
are used to access `data.title` etc. Other variable names, are used from the outer scope.

example: `maxnnz = 1_000; mdlist(@pred(n <= maxnnz))` would produce a list of all
data with less than `maxnnz` structural non-zeros.
"""
macro pred(ex)
    esc(make_pred(ex))
end

