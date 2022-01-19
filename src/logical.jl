
# functions related to the logical operations on patterns
# and predicate functions.
#
"""
    Pattern syntactic sugar

the logical operators `|`, `&`, and `~`
 may be applied to all kind of `Pattern`s with the usual meaning.
+ `~` : unary negation operator - pattern does not match (highest priority)
+ `&` : binary logical and - both patterns match
+ `|` : binary logical or - any of the pattern match (lowest priority)

+ parentheses can be used to overrule operator precedence.

+ `[p...]` is the same as `p[1] | p[2] ...`

+ `(p...)` is the same as `p[1] & p[2] ...`

+ `~(p...) === ~((p...))` - that is `~(p[1] & p[2] ...)`

+ Precedence of '*' is higher that `~` for character and string objects:
 so `~ "a" * "b"  === ~("a" * "b") === ~"ab"`  also `~'a'^2*"b" === ~"aab"`
"""
function logical end

import Base: ~
import LinearAlgebra: isposdef

(&)(p::Tuple, q::Tuple) = tuple(p..., q...)
(&)(p::Tuple, q::Pattern...) = tuple(p..., q...)
(&)(p::Pattern, q::Tuple) = tuple(p, q...)
(&)(p::Pattern, q::Pattern...) = tuple(p, q...)
(|)(p::Pattern, q::Pattern...) = vcat(p, q...)
(~)(c::AbstractChar) = Not(string(c))
(*)(a::Not{<:AbstractString}, b::Union{AbstractString,AbstractChar}) = Not(a.pattern * b)
(~)(p::Not) = p.pattern
Not(p::Not) = p.pattern
(~)(p::Pattern...) = Not(tuple(p...))
(~)(p::Pattern) = p == EMPTY_PATTERN ? ALL_PATTERN : p == ALL_PATTERN ? EMPTY_PATTERN : Not(p)
(~)() = EMPTY_PATTERN
(~)(p::Vector{<:Pattern}) = length(p) == 0 ? ALL_PATTERN : length(p) == 1 ? Not(p[1]) : Not(p)

const EMPTY_PATTERN = Pattern[]
const ALL_PATTERN = ()

###
# simplification of patterns
###
is_pure_string(p::AbstractString) = !occursin(r"[]*?]", p)
is_pure_string(p::Pattern) = false
is_pure_vector(p::Vector) = all(is_pure_string.(p))
is_pure_vector(p::Pattern) = false

###
# Alias resolution
###

function aliasresolve(db::MatrixDatabase, k::AbstractString)
    haskey(db.aliases, k) ? [db.aliases[k]] : String[]
end
function aliasresolve(db::MatrixDatabase, a::Alias{T,<:Integer}) where T
    aliasresolve(db, aliasname(a))
end
function aliasresolve(db::MatrixDatabase, a::Alias{T,<:AbstractVector}) where T
    collect(Iterators.flatten(aliasresolve.(Ref(db), aliasname(a))))
end
aliasresolve(db::MatrixDatabase, a::Alias{RemoteMatrixData{SSRemoteType},Colon}, ) = "*/*"
aliasresolve(db::MatrixDatabase, a::Alias{RemoteMatrixData{MMRemoteType},Colon}) = "*/*/*"
aliasresolve(db::MatrixDatabase, a::Alias{GeneratedMatrixData{:B},Colon}) = :builtin
aliasresolve(db::MatrixDatabase, a::Alias{GeneratedMatrixData{:U},Colon}) = :user

###
# Predefined predicates for MatrixData
###

builtin(p...) = Alias{GeneratedMatrixData{:B}}(p...)
user(p...) = Alias{GeneratedMatrixData{:U}}(p...)
sp1(p...) = Alias{RemoteMatrixData{SSRemoteType}}(p...)
mm1(p...) = Alias{RemoteMatrixData{MMRemoteType}}(p...)
sp(p...) = sp1(p...)
mm(p...) = mm1(p...)

"""
    sp(i, j:k, ...)
    sp(pattern)

The first form with integer and integer range arguments is a pattern selecting
by the id number in the Suite Sparse collection.

The second form is a pattern, which selects a matrix in the Suite Sparse collection, which
corresponds to the pattern by name, even if the name is from the Matrix Market collection.

example:
    mdlist(sp("*/*/1138_bus")) == ["HB/1138_bus"]
"""
sp(p::P) where P<:Pattern = is_ivec(p) ? sp1(p) : Alternate{SSRemoteType,P}(p)

"""
    mm(i, j:k, ...)
    mm(pattern)

The first form with integer and integer range arguments is a pattern selecting
by the id number in the Matrix Market collection.

The second form is a pattern, which selects a matrix in the Matrix Market collection, which
corresponds to the pattern by name, even if the name is from the Suite Sparse collection.

example:
    mdlist(mm("*/1138_bus")) == ["Harwell-Boeing/psadmit/1138_bus"]
"""
mm(p::P) where P<:Pattern = is_ivec(p) ? mm1(p) : Alternate{MMRemoteType,P}(p)

is_ivec(::AbstractVector{<:Integer}) = true
is_ivec(p::AbstractVector) = all(is_ivec.(p))
is_ivec(::typeof(:)) = true
is_ivec(::Any) = false

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
issymmetric(data::MatrixData) = data.name in mdlist(:symmetric)
isskew(data::MatrixData) = false
ishermitian(data::MatrixData) = false
issparse(data::RemoteMatrixData) = true
issparse(data::MatrixData) = data.name in mdlist(:sparse)
isposdef(data::RemoteMatrixData) = hasinfo(data) && data.posdef
isposdef(data::MatrixData) = data.name in mdlist(:posdef)

iscomplex(data::RemoteMatrixData) = _isfield(data, MMFieldComplex)
isreal(data::RemoteMatrixData) = _isfield(data, MMFieldReal)
isinteger(data::RemoteMatrixData) = _isfield(data, MMFieldInteger)
isboolean(data::RemoteMatrixData) = _isfield(data, MMFieldPattern)
iscomplex(data::MatrixData) = false
isreal(data::MatrixData) = false
isinteger(data::MatrixData) = false
isboolean(data::MatrixData) = false

hasinfo(data::RemoteMatrixData) = data.header.m > 0 && data.header.n > 0 # isassigned(data.properties) && data.properties[] !== nothing
hasinfo(data::MatrixData) = false
isremote(data::RemoteMatrixData) = true
isremote(data::MatrixData) = false
isloaded(data::RemoteMatrixData) = hasinfo(data) && !isempty(data.metadata)
isloaded(data::MatrixData) = false
isunloaded(data::RemoteMatrixData) = !isloaded(data)
isunloaded(data::MatrixData) = false
isuser(data::GeneratedMatrixData{:U}) = true
isuser(data::MatrixData) = false
isbuiltin(data::GeneratedMatrixData{:B}) = true
isbuiltin(data::MatrixData) = false
islocal(data::GeneratedMatrixData) = true
islocal(data::MatrixData) = false

"""
    pred(f::Function, s::Symbol...)

Return a predicate function, which assigns to a each `data::MatrixData`
iff
*    `hasdata(data)` and
*    all symbols `s` are property names of `data` and
*    `f` applied to the tuple of values of those properties returns `true`
"""
function pred(f::Function, s::Symbol...)
    data::MatrixData -> hasinfo(data) &&
    s ⊆ propertynames(data) &&
    f(getproperty.(Ref(data), s)...)
end

"""
    prednzdev(deviation)

Test predicate - does number of stored (structural) non-zeros deviate from nnz
by more than `deviation`. That would indicate a data error or high number of stored zeros. 
"""
function prednzdev(dev::AbstractFloat=0.1)
    function f(data::RemoteMatrixData)
        n1, n2 = extremnnz(data)
        n1 -= Int(floor(n1 * dev))
        n2 += Int(floor(n2 * dev))
        isloaded(data) && ! ( n1 <= data.nnz <= n2 )
    end
    f(::MatrixData) = false
    f
end

"""
    keyword(word::Union{AbstractString,Tuple,Vector})

Predicate function checks, if `word` is contained in on to the textual
metadata fields `[:notes, :title, :kind, :author]`.
Tuples and Vectors are interpreted as `AND` resp. `OR`.
"""
function keyword(s::AbstractString)
    function f(data::RemoteMatrixData)
        text = join([data.notes, data.title, data.author, data.kind], ' ')
        match(Regex("\\b$s\\b", "i"), text) !== nothing
    end
    f(::MatrixData) = false
    f
end
# this is to translate `keyword("a" & "b")` to `keyword("a") & keyword("b")`
keyword(t::Tuple) = Tuple(keyword.(t))
keyword(t::AbstractVector{<:AbstractString}) = [keyword(x) for x in t]

"""
    hasdata(meta::Union{Symbol,Tuple,Vector})

Predicate function checks, if matrix data have metadata symbol `meta`.
Tuples and Vectors are interpreted as `AND` resp. `OR`.
"""
function hasdata(s::Symbol, s2::Symbol...)
    function f(data::MatrixData)
        ms = metasymbols(data)
        s in ms && issubset(s2, metasymbols(data))
    end
    f
end
hasdata(t::NTuple{N,Symbol} where N) = hasdata(t...)
hasdata(t::Tuple) = Tuple(hasdata.(t))
hasdata(t::AbstractVector) = [hasdata(x) for x in t]

"""
    check_symbols(p::Pattern)

throw `ArgumentError` if pattern uses unknown symbol as a group name.
"""
function check_symbols(p::Pattern)
    s = setdiff(filter(x->x isa Symbol, flatten_pattern(p)), listgroups())
    isempty(s) || argerr("The following symbols are no group names: $s")
end

############################
# Predicate generating macro
############################
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

"""
    PROPS

The symbols, which are naming properties of `MatrixData`
"""
const PROPS = (:name, :id, :metadata, fieldnames(MetaInfo)...)

function make_pred(ex)
    syms = extract_symbols(ex) ∩ PROPS 
    :(MatrixDepot.pred($(make_func(syms, ex)), $(QuoteNode.(syms)...)))
end

"""
    @pred(expression)

Generate a predicate function using the expression as function body. Variable names
within the expression, which are properties of `RemoteMatrixData` (e.g. `title`, `m`, `nnz`)
are used to access `data.title` etc. Other variable names, are used from the outer scope.

example: `maxnnz = 1_000; listnames(@pred(n <= maxnnz))` would produce a list of all
data with less than `maxnnz` structural non-zeros.
"""
macro pred(ex)
    esc(make_pred(ex))
end

