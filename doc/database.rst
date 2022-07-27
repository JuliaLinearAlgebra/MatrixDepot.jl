Basics
=========

MatrixDepot is a collection of matrix centric problems from different sources.
Locally defined problems are generated in the local machine. They need parameters, for example the dimensions, to specify a matrix.
Remote problems are accessing one of 2 remote repositories, the Suite Sparse Matrix Collection, and the Matrix Market Collection.
In both cases, the matrices and additional data are cached on the local disk.

Matrix Names
------------

Each problem has a canonical name in the form of a relative pathname.
For the local problems, the name contains no slash `/`. Example: `"vand"`.
For the Suite Sparse Collection, the name contains one slash. Example: `"HB/ash219"`.
For the Matrix Market Collection, the name contains two slashes. Example: `"Harwell-Boeing/smtape/ash219"`.

Each canonical name is mapped to a `MatrixData` object by the call `data = MatrixDepot.mdata(name::String)`.
Reversely, the name can be recovered by `name = MatrixDepot.filename(data::MatrixData)`.

Besides the name, within each category, the problems have integer identifiers, counting from 1. (`data.id`)

Patterns
--------

A pattern is a way to define a subset of the matrices contained in the database.
While patterns were originally restricted to search patterns in the matrix names, they allow now filtering by many matrix properties.

Patterns are typically used as arguments to `mdlist(pattern)`, which returns a vector of known matrix names matching the pattern.
Each pattern can be mapped to the characteristic function of the subset, it represents via `f = charfun(pattern)`. This way `f(name)` returns
true iff `name` is element of the subset.

Name patterns:

- Pure strings (not containing special characters `'*'`, `'?'`, `'['`, `']'`),
  exactly match a canonical name.

- Shell pattern strings (the special characters are interpreted like in a shell, `"**"` also matches `'/'`)

- Regular expressions (e.g. `r"^v.*d$"`, applied to the name)

- Group names of all defined groups `listgroups()`, as a symbol

- Integer filters (`builtin(1)`, `user(:)`, `sp(2800:2810)`, `mm(1,2,5:7)`), using integer identifiers, ranges, colon

- Forcing collection filters (`sp("H*/*/1*")`, `mm("HB/1*")`), which find the corresponding matrix of the alternative collection, if that exists.

- Predefined predicate functions (actually every function `f(::MatrixData)::Boolean` works.)

- -  isboolean, isinteger, isreal, iscomplex -- (element type - boolean == pattern)

- -  isgeneral, issymmetric, ishermitian, isskew -- (matrix symmetry)

- -  isbuiltin, isuser, islocal, isremote, isloaded, isunloaded -- (status)

- -  issvdok -- (are svd data available?)

- -  keyword(string expression) -- (single keyword or keywords connected with `&` or `|`)

- -  hasdata(symbol) -- (are data with that name available? `mdlist(hasdata(:b))` --> `mdopen(...).b`)

- -  @pred(expression) -- (expressions can access all matrix properties, example `@pred(1 < m < 90)`)

@pred Function
--------------

The `@pred` macro defines a predicate function by an expression, which is applied to the `data` object.
All names, which are found in `propertynames(data)` are referring to those properties of `data`.
Other names refer to global names, which are accessible from the definition of the `@pred` function.

Logical Expressions
-------------------

Patterns of all kinds may be combined using the logical operators `~`, `&`, `|`.
While that may be considered as syntactical sugar (for `Not()`, `(,)`, `[,]`) it is a case of type piracy,
because `Pattern` includes the base types of `String`, `RegEx`, `Tuple`, `Vector`, and `Function`.

In a future release, these operators will need to be opted in by `using MatrixDepot.Logicals`.

Examples::

    julia> data = MatrixDepot.mdata("Sorensen/Linux_call_graph")
    (IG Sorensen/Linux_call_graph(#2654)  324085x324085(1208908) 2013 [A] 'directed weighted graph' [Call graph of the Linux 3.7.10 kernel])

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
