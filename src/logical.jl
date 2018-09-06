
# functions related to the logical operations on patterns
# and predicate functions.
#
export builtin, user, uf, tamu, mm, ¬, logical

export isgeneral, issymmetric, isskew, ishermitian
export iscomplex, isreal, isinteger, ispattern
export isremote, islocal, isloaded, isunloaded, isbuiltin, isuser
export predm, predn, prednz, predmn

import Base: isreal, isinteger
import LinearAlgebra: issymmetric, ishermitian

import Base: &, |, *

"""
    Syntactic sugar

the logical operators `|`, `&`, and `¬` (`\\neg`)
 may be applied to all kind of `Pattern`s with the usual meaning.
+ `¬` : unary negation operator - pattern does not match (highest priority) 
+ `&` : binary logical and - both patterns match
+ `|` : binary logical or - any of the pattern match (lowest priority)

+ The operator `*` for string concatenation is treated special
 so `¬ "a" * "b"  === ¬("a" * "b") === ¬"ab"`  also `¬'a'^2*"b" === ¬"aab"`
+ also `¬(p1...) === ¬((p...))` has been implemented.
"""
logical() = "¬|?"

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

function aliasresolve(k::AbstractString, db::MatrixDatabase)
    haskey(db.aliases, k) ? [db.aliases[k]] : String[]
end
function aliasresolve(a::Alias{T,<:Integer}, db::MatrixDatabase) where T
    aliasresolve(aliasname(a), db)
end
function aliasresolve(a::Alias{T,<:AbstractVector{<:Integer}}, db::MatrixDatabase) where T
    aliasr2(x) = aliasresolve(x, db)
    flatten(aliasr2.(aliasname(a)))
end
aliasresolve(a::Alias{RemoteMatrixData{TURemoteType},Colon}, db::MatrixDatabase) = "*/*"
aliasresolve(a::Alias{RemoteMatrixData{MMRemoteType},Colon}, db::MatrixDatabase) = "*/*/*"
aliasresolve(a::Alias{GeneratedBuiltinMatrixData,Colon}, db::MatrixDatabase) = :builtin
aliasresolve(a::Alias{GeneratedUserMatrixData,Colon}, db::MatrixDatabase) = :user

###
# Predefined predicates for MatrixData
###

builtin(p::IntOrVec) = Alias{GeneratedBuiltinMatrixData}(p)
user(p::IntOrVec) = Alias{GeneratedUserMatrixData}(p)
uf(p::IntOrVec) = Alias{RemoteMatrixData{TURemoteType}}(p)
tamu(p::IntOrVec) = Alias{RemoteMatrixData{TURemoteType}}(p)
mm(p::IntOrVec) = Alias{RemoteMatrixData{MMRemoteType}}(p)

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

hasinfo(data::RemoteMatrixData) = data.header.m > 0 && data.header.n > 0
hasinfo(data::MatrixData) = false
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

predm(f::Function) = data::MatrixData -> hasinfo(data) && f(row_num(data))
predn(f::Function) = data::MatrixData -> hasinfo(data) && f(col_num(data))
prednz(f::Function) = data::MatrixData -> hasinfo(data) && f(nz_num(data))
preddnz(f::Function) = data::MatrixData -> hasinfo(data) && f(dnz_num(data))
function predmn(f::Function)
    data::MatrixData -> hasinfo(data) && f(row_num(data), col_num(data))
end
