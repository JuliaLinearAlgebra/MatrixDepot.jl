
# ![logo](doc/logo2.png) Matrix Depot

[![Build Status][gha-img]][gha-url]    [![Coverage Status][codecov-img]][codecov-url]    [![Documentation Status][rtd-img]][rtd-url]

An extensible test matrix collection for Julia.

* [Release Notes](https://github.com/JuliaMatrices/MatrixDepot.jl/blob/master/NEWS.md)


Give access to a wealth of sample and test matrices and accompanying data.
A set of matrices is generated locally (with arguments controlling the special case).
Another set is loaded from one of the publicly accessible matrix collections
[`SuiteSparse Matrix Collection`](https://sparse.tamu.edu) (formerly `University of Florida Matrix Collection`)
and the [`Matrix Market Collection`](https://math.nist.gov/MatrixMarket), the latter being obsolescent.

Access is like

```julia
using MatrixDepot
?MatrixDepot                    # display package help info
A = matrixdepot("hilb", 10)     # locally generated hilbert matrix dimensions (10,10)
A = matrixdepot("HB/1138_bus")  # named matrix of the SuiteSparse Collection
```

or

```julia
md = mdopen("*/bfly")   # named matrix with some extra data
A = md.A
co = md.coord
tx = md("Gname_10.txt")
md.<tab><tab>           # overview of the "fields" of md returning like
                        # A m n dnz nnz coord Gname_10.txt  G_10 Gcoord_10
```

or also

```julia
mdinfo("gravity")                 # text info about the selected matrix
md = mdopen("gravity", 10, false) # locally generated example with rhs and solution
A = md.A
b = md.b
x = md.x
```

## Install

**NOTE:** If you use Windows, you need to install MinGW/MSYS or Cygwin
in order to use the SuiteSparse sparse and MatrixMarket matrix collection interface.

To install the release version, type

```julia
julia> Pkg.add("MatrixDepot")
```

## Usage

### Naming

#### Matrix Names

Every Matrix type has a unique name, which is a string of one of the forms:

  1. `"name"` - used for matrices, which are generated locally.
  2. `"dir/name"` - for all matrices of the `SuiteSparse` collection.
  3. `"dir/subdir/name"` - for all matrices of the `MatrixMarket` collection.

The names are similar to relative path names, separated by a slash character.
The components of the name must not contain any of the characters `"/*[]"`.

#### Groups

a set of matrices may be assigned to predefined or user-defined groups.
The group names are represented as `Julia` symbols in the form `:symmetric`.
The group names are therefore restricted to valid `Julia` identifiers, that means
start with a letter and contain only letters, digits, and `'_'`.

#### Numeric Identifiers

Every matrix has a numeric identifier, which is unique for its area:

* `builtin(id)` - one of the built-in matrix generators - currently `id ∈ 1:59`.

* `user(id)` - a user-defined matrix generator - starting with `1`.

* `sp(id)` - one of the `SuiteSparse` collection. The integer ids are the
  'official' ident numbers assigned by the collection. Currently `id ∈ 1:3000`.

* `mm(id)` - one of the `MatrixMarket` collection. Here id follows the ordering
  of the index file of the collection.

### Sets of Matrix Names - Pattern

For some functions it makes sense to have lists of matrix names to operate on, for
example to select a set matrices with certain properties. These sets are described
by 'Patterns', which are applied to matrix names and also to other matrix properties.

The following pattern types are supported:

  1. `"name"` - a string matching exactly a matrix name
  2. `"shell-pattern"` - a string with shell wildcards `'?', '*', "[...]"` included.

  3. `r"egular expression"` - a regular expression to match the matrix name.

  4. `:group` - one of the defined group names; match all matrices in the group

  5. `qualified numeric identifiers` - examples `builtin(10)`, `sp(1:5, 7)`, `mm(1), sp(:)`

  6. `predicate_function` - the name of a predefined or user-defined boolean function of the internal data type `MatrixData`. Example: `issymmetric`.

  7. `abstract vector of sub-patterns` - `OR` - any of the sub-pattern matches

  8. `tuple of sub-patterns` - `AND` - all of the sub-patterns match

  9. `~pattern` - negation of a pattern the \neg - operator ~ may be applied to all patterns

To express `OR` and `AND`, the binary operators `|` and `&` and `(` / `)` are preferred.

Examples:

```julia
"gravity" | "HB/*" & ~(ishermitian & iscomplex) & ~sp(20:30)
```

The set of all known matrices can be expressed as empty tuple `()`. In a shell-
pattern the double `**` matches also slash characters, in contrast to the single `*`.

A convenient form of a predicate-generator is

```julia
@pred(expression)
```

where expression is a valid `Julia` boolean expression, which may access all
properties of `MatrixData` as literal variable names.

Examples:

`@pred(author == "J. Brown")` is translated to:
`d -> :author in propertynames(d) && d.author == "J. Brown"`

`@pred(500_000 <= n * m < 1_000_000)` restricts the size of matched matrices.

`@pred(10^4 <= n <= 2*10^4 && n == m && nnz / n > 10 )` in average more than 10 entries per row

There is s set of predefined predicate functions including:
`(issymmetric, ishermitian, isgeneral, isskew, isreal, iscomplex, isboolean,
islocal, isremote, isloaded, isunloaded, isbuiltin, isuser, issparse)`

Special predicate generators `keyword(word...)` and `hasdata(symbol...)` allow to
support keyword-search and check for the existence of meta-data.
For example: `hasdata(:x) & ~keyword("fluid"` provides solution (x) and does not mention "fluid".

### Number of nonzeros

Beware that some sparse matrices contain non-structural zeros, that is, coefficients stored explicitly but whose value is `0`.
In this case a discrepancy between nnz(A) and sum(!iszero, A) will be observed.

## Accessing Data

### Listing matrices with certain properties

```julia
mdinfo()           # overview
listgroups()       # list all defined group names
mdlist(pattern)    # array of matrix names according to pattern
listdata(pattern)  # array of `MatrixData`objects according to pattern
listnames(pattern) # MD-formatted listing of all names according to pattern
listdir("*//*") # MD-formatted -  group over part before `//` - count matching
```

### Information about matrices

```julia
mdinfo()        # overview over database
mdinfo(pattern) # individual documentation about matrix(es) matching pattern
```

### Generating a matrix

`A = matrixdepot("kahan", 10)` generates a matrix using one of the built-in generators

`md = mdopen("kahan", 10)` returns a handle `md`; matrix can be obtained by
`A = md.A`

### Accessing Meta-Data

In general the first form is preferable, if only the pure matrix is required.
For remote collections no arguments are used.

The second form allows to access all types of 'meta-data', which may be available for some local or remote matrices.

Examples:

`md = mdopen("spikes", 5, false); A = md.A; b = md.b; x = md.x`

`md = mdopen("Rommes/bips07_1998"); A = md.A; v = md.iv; title = md.data.title;
 nodenames = md("nodename.txt")`

The last example shows, how to access textual meta-data, when the name contains
`Julia` non-word characters. Also if the metadata-name is stored in a variable,
the last form has to be used.

`meta = metasymbols(md)[2]; sec_matrix = md(meta)`

The function `metasymbols` returns a list of all symbols denoting metadata
provided by `md`. Whether expressed as symbols or strings does not matter.

The system function `propertynames(md)` returns all data of `md`. That includes
size information and metadata.

`propertynames(md.data)` gives an overview about all attributes of the
`md.data::MatrixData`, which can for example be used in the `@pred` definitions.

### Backoffice Jobs

The remote data are originally stored at the remote web-site of one of the
matrix collections. Before they are presented to the user, they are downloaded
to local disk storage, which serves as a permanent cache.

By default, the data directory is a scratchspace managed by [`Scratch.jl`](https://github.com/JuliaPackaging/Scratch.jl), but can be changed by setting the `MATRIXDEPOT_DATA` environment variable.

The data directory can be queried by

    julia> MatrixDepot.data_dir()
    "/home/.../.julia/scratchspaces/b51810bb-c9f3-55da-ae3c-350fc1fbce05/data

The occasional user needs not bother about downloads, because that is done in
the background if matrix files are missing on the local disk.

The same is true for the data required by `mdinfo(pattern)`. Actually these are
stored in separate files if the full matrix files (which may be huge) are not yet loaded.

#### Bulk Downloads

##### Load Header Data

A download job to transmit a subset of remote matrix files may be started to
load header data for all files. Header data always include the matrix type
according to the matrix-market-format and the size values `m` row-number,
`n` = columns-number, and `dnz` number of stored data of the main sparse matrix.

`MatrixDepot.loadinfo(pattern)` where `pattern` defines the subset.

That is possible for the [SuiteSparse collection](https://sparse.tamu.edu) and the
[NIST MatrixMarket collection](https://math.nist.gov/MatrixMarket).
The patterns can always refer to matrix names and id numbers.
In the case of `SuiteSparse` collection, also the metadata
`"date"`, `"kind"`, `"m"`, `"n"`, `"nnz"`
are available and can be used, before individual matrix data
have been loaded. They are contained in a data file obtained from the remote site.
For `MatrixMarket` collection, patterns are restricted to names and id numbers.

In general it would be possible by `loadinfo("**")` to load all header data. That
would last maybe an hour and generate some traffic for the remote sites.
Nevertheless it is not necessary to do so, if you don't need the header data
for the following task.

##### Load Complete Data Files

**`MatrixDepot.load(pattern)`** loads all data files for the patterns.
Patterns can only refer to attributes, which are already available.
In the case of `SuiteSparse` that includes the size info `"date"`, `"kind"`,
`"m"`, `"n"`, and `"nnz"` and all additional attributes loaded in the previous step,
which include `"author"`, `"title"`, `"notes"`, and keywords.
In the case of `MatrixMarket` you can only refer to `"m"`, `"n"`, and `"dnz"`,
if previously loaded with the header data.

Please do not:
`MatrixDepot.load("**")`. That would require some day(s) to finish and include
some really big data files (~100GB), which could be more than your disks can hold.

Make a reasonable selection, before you start a bulk download.
Local and already loaded matrices are skipped automatically.

Example:

`MatrixDepot.load(sp(:) & @pred(nnz < 100_000))` to download only problems with given
number of stored entries in the main matrix.

## Sample Session

To see an overview of the matrices in the collection, type

```julia
julia> using MatrixDepot

julia> mdinfo()
  Currently loaded Matrices
  –––––––––––––––––––––––––––

builtin(#)
–––––––––– ––––––––––– ––––––––––– ––––––––––– –––––––––– –––––––––––– ––––––––––– ––––––––––– ––––––––––––– ––––––––––––
1 baart    7 circul    13 fiedler  19 gravity  25 invhilb 31 magic     37 parter   43 randcorr 49 shaw       55 ursell
2 binomial 8 clement   14 forsythe 20 grcar    26 invol   32 minij     38 pascal   44 rando    50 smallworld 56 vand
3 blur     9 companion 15 foxgood  21 hadamard 27 kahan   33 moler     39 pei      45 randsvd  51 spikes     57 wathen
4 cauchy   10 deriv2   16 frank    22 hankel   28 kms     34 neumann   40 phillips 46 rohess   52 toeplitz   58 wilkinson
5 chebspec 11 dingdong 17 gilbert  23 heat     29 lehmer  35 oscillate 41 poisson  47 rosser   53 tridiag    59 wing
6 chow     12 erdrey   18 golub    24 hilb     30 lotkin  36 parallax  42 prolate  48 sampling 54 triw

user(#)
–––––––––
1 randsym

Groups
–––––– ––––––– ––––– –––– ––––– ––––– ––––––– ––––––– –––––– –––––– ––––––– –––––– –––––––––
all    builtin local user eigen graph illcond inverse posdef random regprob sparse symmetric

Suite Sparse of
–––––––––––– ––––
2770         2833

MatrixMarket of
–––––––––––– –––
488          498

```

We can generate a 4-by-4 Hilbert matrix by typing

```julia
julia> matrixdepot("hilb", 4)
4x4 Array{Float64,2}:
 1.0       0.5       0.333333  0.25
 0.5       0.333333  0.25      0.2
 0.333333  0.25      0.2       0.166667
 0.25      0.2       0.166667  0.142857
```

We can type the matrix name to get documentation about the matrix.

```julia
julia> mdinfo("hilb")
     Hilbert matrix
    ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

  The Hilbert matrix has (i,j) element 1/(i+j-1). It is notorious for being
  ill conditioned. It is symmetric positive definite and totally positive.

  Input options:

    •  [type,] dim: the dimension of the matrix;

    •  [type,] row_dim, col_dim: the row and column dimensions.

  Groups: ["inverse", "ill-cond", "symmetric", "pos-def"]

  References:

  M. D. Choi, Tricks or treats with the Hilbert matrix, Amer. Math. Monthly,
  90 (1983), pp. 301-312.

  N. J. Higham, Accuracy and Stability of Numerical Algorithms, second
  edition, Society for Industrial and Applied Mathematics, Philadelphia, PA,
  USA, 2002; sec. 28.1.
```

We can also specify the data type for locally generated matrices.

```julia
julia> matrixdepot("hilb", Float16, 5, 3)
5x3 Array{Float16,2}:
 1.0      0.5      0.33325
 0.5      0.33325  0.25
 0.33325  0.25     0.19995
 0.25     0.19995  0.16663
 0.19995  0.16663  0.14282

julia> matrixdepot("hilb", Rational{Int}, 4)
4x4 Array{Rational{T<:Integer},2}:
 1//1  1//2  1//3  1//4
 1//2  1//3  1//4  1//5
 1//3  1//4  1//5  1//6
 1//4  1//5  1//6  1//7
```

Matrices can be accessed by a variety of patterns and composed patterns.
Integer numbers `i` refer to the ident numbers in `sp(i)`, `mm(i)`, `builtin(i)`, `user(i)`.
Here `sp` ... denote the supported matrix collections SuiteSparse (formerly UFL),
Matrix Market, built-in, user-defined.

```julia
julia> mdlist(sp(1))    # here sp(1) is the ident number of the SuiteSparse collection
list(1)
–––––––––––
HB/1138_bus

julia> listnames(builtin(1, 5:10))    # the internal numbering of the builtin-functions
list(7)
––––––– –––––––– –––– –––––– ––––––– ––––––––– ––––––
baart   chebspec chow circul clement companion deriv2

julia> mdlist(builtin(1:4, 6, 10:15) | user(1:10) )
12-element Array{String,1}:
 "baart"
 "binomial"
 "blur"
 "cauchy"
 "chow"
 "deriv2"
 "dingdong"
 "erdrey"
 "fiedler"
 "forsythe"
 "foxgood"
 "randsym"
```

While the `listnames` command renders the output as markdown table, the internal
`mdlist` produces an array of valid matrix names.

We can type a group name to see all the matrices in that group. Group names are
always written as symbols to distinguish them form matrix names and pattern, which
are always strings.

```julia
julia> listnames(:symmetric)
list(22)
–––––––– –––––––– ––––––– –––––– ––––––––– –––––––– ––––––– –––––––––
cauchy   dingdong hilb    lehmer oscillate poisson  randsym wilkinson
circul   fiedler  invhilb minij  pascal    prolate  tridiag
clement  hankel   kms     moler  pei       randcorr wathen
```

## Extend Matrix Depot

It is possible to extend the builtin local problems with user defined generators and groups.
We can add [new matrix generators](http://matrix-depot.readthedocs.org/en/latest/user.html)
and define [new groups of matrices](http://matrix-depot.readthedocs.org/en/latest/properties.html).

## References

* Weijian Zhang and Nicholas J. Higham,
  "Matrix Depot: An Extensible Test Matrix Collection for Julia",
  *PeerJ Comput. Sci.*, 2:e58 (2016),
  [[pdf]](https://peerj.com/articles/cs-58/)

* Nicholas J. Higham,
  "Algorithm 694, A Collection of Test Matrices in MATLAB",
  *ACM Trans. Math. Software*,
  vol. 17. (1991), pp 289-305
  [[pdf]](http://www.maths.manchester.ac.uk/~higham/narep/narep172.pdf)
  [[doi]](https://dx.doi.org/10.1145/114697.116805)
* T.A. Davis and Y. Hu,
  "The University of Florida Sparse Matrix Collection",
  *ACM Transaction on Mathematical Software*,
  vol. 38, Issue 1, (2011), pp 1:1-1:25
  [[pdf]](http://www.cise.ufl.edu/research/sparse/techreports/matrices.pdf)

* R.F. Boisvert, R. Pozo, K. A. Remington, R. F. Barrett, & J. Dongarra,
  " Matrix Market: a web resource for test matrix collections",
  *Quality of Numerical Software* (1996) (pp. 125-137).
  [[pdf]](http://www.netlib.org/utk/people/JackDongarra/pdf/matrixmarket.pdf)

* Per Christian Hansen,
  "Test Matrices for Regularization Methods",
  *SIAM Journal on Scientific Computing*,
  vol. 16, 2, (1995) pp.506-512.
  [[pdf]](http://epubs.siam.org/doi/abs/10.1137/0916032)
  [[doi]](https://dx.doi.org/10.1137/0916032)

[gha-img]: https://github.com/JuliaMatrices/MatrixDepot.jl/workflows/CI/badge.svg
[gha-url]: https://github.com/JuliaMatrices/MatrixDepot.jl/actions?query=workflow%3ACI

[codecov-img]: https://codecov.io/gh/JuliaMatrices/MatrixDepot.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaMatrices/MatrixDepot.jl

[rtd-img]: https://readthedocs.org/projects/matrix-depot/badge/?version=latest
[rtd-url]: https://matrix-depot.readthedocs.io/en/latest/?badge=latest
