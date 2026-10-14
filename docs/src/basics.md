# Basics

```@setup basics
using MatrixDepot
```

MatrixDepot is a collection of matrix-centric problems from different sources.
*Local* problems are generated on the local machine and take parameters, for example
the dimension. *Remote* problems come from one of two repositories, the SuiteSparse
Matrix Collection and the Matrix Market Collection. Downloaded matrices and additional
data are cached on the local disk.

The internal database is loaded automatically by `using MatrixDepot`, which downloads the
index files of the remote collections on first use.

## Matrix names

Each problem has a canonical name in the form of a relative path name:

| collection    | form                  | example                          |
|:--------------|:----------------------|:---------------------------------|
| local         | `"name"`              | `"vand"`                         |
| SuiteSparse   | `"dir/name"`          | `"HB/ash219"`                    |
| Matrix Market | `"dir/subdir/name"`   | `"Harwell-Boeing/smtape/ash219"` |

Name components must not contain any of the characters `/*[]`. Within each collection,
problems also have integer identifiers counting from 1, accessible as
[`builtin`](@ref), [`user`](@ref), [`sp`](@ref) and [`mm`](@ref).

Each canonical name is mapped to a `MatrixData` object by `MatrixDepot.mdata(name)`.

## Generating matrices

[`matrixdepot`](@ref) returns the matrix selected by a name or a pattern matching exactly
one problem; further arguments are passed to the generator:

```@example basics
matrixdepot("hilb", 5, 4)
```

An element type can be given as the first parameter:

```@example basics
matrixdepot("hilb", Float16, 3)
```

[`mdinfo`](@ref) shows the documentation of a problem, including its input options,
groups and references:

```@example basics
mdinfo("moler")
```

[`mdopen`](@ref) returns a descriptor that gives access to the matrix and accompanying data.
For the regularization problems this includes a right-hand side and solution:

```@example basics
md = mdopen("gravity", 4, false)
md.b
```

Called without arguments, `mdinfo()` gives an overview of all loaded problems, groups
and remote collections.

## Patterns

A pattern selects a subset of the problems in the database. Patterns are accepted by
[`matrixdepot`](@ref), [`mdinfo`](@ref), [`mdlist`](@ref), [`listnames`](@ref),
[`listdata`](@ref), [`mdopen`](@ref) (single match required), `MatrixDepot.load`,
`MatrixDepot.loadinfo` and `MatrixDepot.loadsvd`. [`listdir`](@ref) accepts only strings.
[`mdlist`](@ref) returns the names matching a pattern, and [`charfun`](@ref) turns a
pattern into the characteristic function of the subset.

The elementary patterns are:

1. A name, or a shell pattern with the wildcards `?` (one character except `/`),
   `*` (any sequence except `/`), `**` (any sequence including `/`) and `[...]`.
   Example: `"*/???/w*"`, `"**/1138*"`.
2. A regular expression applied to the name, e.g. `r"^v.*d$"`.
3. An integer identifier `builtin(n)`, `user(n)`, `sp(n)` or `mm(n)`, where `n` may be
   integers, ranges or `:`. Example: `builtin(1, 3, 10:11)`.
4. An alternate-collection pattern `sp(pattern)` or `mm(pattern)`, which selects the
   corresponding matrix of the other remote collection, if it exists (names are
   sometimes mangled). Example: `sp("*/*/1138_bus")` selects `"HB/1138_bus"`.
5. A group name as a `Symbol`, e.g. `:symmetric`; see [Groups](groups.md).
6. A predicate function `f(::MatrixData)::Bool`, for example
   `data::MatrixDepot.RemoteMatrixData -> data.n <= 100`, or the shorthand [`@pred`](@ref):
   `@pred(n <= 100)`.
7. A predefined predicate:
   * element type: `isboolean`, `isinteger`, `isreal`, `iscomplex`
   * symmetry: `isgeneral`, `issymmetric`, `ishermitian`, `isskew`
   * status: `isbuiltin`, `isuser`, `islocal`, `isremote`, `isloaded`, `isunloaded`
   * properties: `issvdok` (SVD data available), `isposdef`
   * [`keyword`](@ref)`(expr)`: single keyword or keywords combined with `&` and `|`
   * [`hasdata`](@ref)`(symbol)`: data with that name are available, e.g.
     `mdlist(hasdata(:b))` lists problems with `mdopen(...).b`.

The [`@pred`](@ref) macro defines a predicate by an expression applied to the `data`
object. Names found in `propertynames(data)` refer to those properties; other names refer
to global names visible where `@pred` is used.

### Combining patterns

Patterns of all kinds combine with the operators `|` (or), `&` (and) and `~` (not);
equivalently, a vector of patterns means *or* and a tuple means *and*.
See also [`logical`](@ref).

```@example basics
mdlist((issymmetric | ishermitian) & :illcond & ~:posdef)
```

```@example basics
listnames(:symmetric & :inverse & :illcond & :posdef)
```

```@example basics
mdlist(builtin(1:4, 6, 10:15))
```

Patterns over remote properties work the same way:

```julia-repl
julia> mdlist(@pred(324085 <= n < 325700))
2-element Vector{String}:
 "LAW/cnr-2000"
 "Sorensen/Linux_call_graph"

julia> mdlist("*/*/*" & @pred(n > 50000))
2-element Vector{String}:
 "misc/cylshell/s3dkq4m2"
 "misc/cylshell/s3dkt3m2"

julia> mdlist(sp(:)) == mdlist("*/*")
true
```
