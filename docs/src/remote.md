# Remote Collections

Remote matrices are downloaded on first access and cached on the local disk.
[`matrixdepot`](@ref) and [`mdopen`](@ref) download the full matrix data implicitly,
[`mdinfo`](@ref) downloads only the header data. A bulk of data can be fetched in advance
with `MatrixDepot.loadinfo(pattern)` and `MatrixDepot.load(pattern)`; SVD data have to be
requested explicitly with `MatrixDepot.loadsvd(pattern)`.

## SuiteSparse Matrix Collection

Problems of the [SuiteSparse Matrix Collection](https://sparse.tamu.edu) are named
`"collection_name/matrix_name"`:

```julia-repl
julia> md = mdopen("SNAP/web-Google")
(PG SNAP/web-Google(#2301)  916428x916428(5105039) 2002 [A] 'Directed Graph' [Web graph from Google]()

julia> M = md.A
916428×916428 SparseMatrixCSC{Bool, Int64} with 5105039 stored entries:
...
```

If a matrix name is unique, `matrixdepot("**/matrix_name")` works as well; an error is
thrown if the pattern matches more than one matrix. `listnames("*/*")` lists the whole
collection.

The available data and metadata are shown by the descriptor and by [`metasymbols`](@ref):

```julia-repl
julia> md = mdopen("TKK/t520")
(IS TKK/t520(#1908)  5563x5563(286341/145952) 2008 [A, b, coord] 'Structural Problem' [T-beam, L = 520 mm, Quadratic four node DK type elements.  R Kouhia]()

julia> metasymbols(md)
(:A, :b, :coord)

julia> md.coord
5563×3 Matrix{Float64}:
 ...
```

Further properties of the problem, such as `rank`, `cond` or the singular values `sv`,
are fields of `md.data`. The predicate `issvdok` selects problems with SVD data:

```julia-repl
julia> mdlist(issvdok & @pred(5900 <= n < 6000))
8-element Vector{String}:
 "AG-Monien/ukerbe1"
 "Cote/mplate"
 "HB/man_5976"
 "Hamrle/Hamrle2"
 "JGD_Homology/cis-n4c6-b4"
 "JGD_Homology/n4c6-b4"
 "Schenk_IBMNA/c-32"
 "TOKAMAK/utm5940"

julia> md = mdopen("AG-Monien/ukerbe1");

julia> md.data.rank
4108
```

For the meaning of the property names see <https://sparse.tamu.edu/statistics>.

## NIST Matrix Market

Problems of the [Matrix Market](https://math.nist.gov/MatrixMarket/) are named
`"collection_name/set_name/matrix_name"`, for example `"Harwell-Boeing/lanpro/nos5"`.
As always, components may be replaced by wildcards as long as the pattern matches one
matrix:

```julia-repl
julia> md = mdopen("*/*/bp__1400")
(RG Harwell-Boeing/smtape/bp__1400(#M93)  822x822(4790)  [A] '' []()

julia> mdinfo(md)
  Harwell-Boeing/smtape/bp__1400
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

  MatrixMarket matrix coordinate real general

  822 822 4790
```

There is typically no header information in this collection besides `m`, `n` and `dnz`.
