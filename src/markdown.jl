
using Markdown

"return information about all matrices selected by pattern"
mdinfo(p::Pattern) = mdinfo(MATRIX_DB, p)
function mdinfo(db::MatrixDatabase, p::Pattern)
    check_symbols(p)
    mdbuffer = Markdown.MD([])
    md = mdbuffer
    for name in mdlist(p)
        try
            md = mdinfo(db.data[name])
        catch
            md = Markdown.parse("# $name\nno info available")
        finally
            append!(mdbuffer.content, md.content)
        end
    end
    mdbuffer
end

mdinfo(md::MatrixDescriptor) = mdinfo(md.data)

_mdheader(md::Markdown.MD, p, o) = isempty(md.content) ? (nothing, o) : _mdheader(md.content[1], md, o)
_mdheader(md::Markdown.Header, p, o) = (md, p)
_mdheader(md, p, o) = (nothing, o)

function mdinfo(data::GeneratedMatrixData)
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

function mdinfo(data::RemoteMatrixData)
    file = verify_loadinfo(data)
    txt = mmreadcomment(file)
    txt = _repl(txt, r"^%-+$"m => "---", r"^%%" => "###### ", r"%+" => "* ")
    md = Markdown.parse(txt)
    insert!(md.content, 1, Markdown.Header{1}([data.name]))
    md
end

"""
format output list to multi- column table form adapted to display width.
"""
function buildnametable(hdr, name::AbstractVector, maxrow::Integer=0)
    maxrow = maxrow <= 0 ? displaysize(stdout)[2] + maxrow : maxrow
    maxrow = clamp(maxrow, 20, 1000)
    hdr isa Vector || ( hdr = [hdr] )
    n = length(name)
    name = string.(name)
    cols = 1
    while sumwidth(name, cols+1) <= maxrow - 1 * cols && cols < n
        cols += 1
    end
    rows = (length(name) + cols - 1) ÷ cols
    n = length(name)
    while n < rows * cols
        push!(name, "")
        n += 1
    end
    name = reshape(name, rows, cols)
    name = vecvec(name)
    insert!(name, 1, [string.(hdr)..., tab("", cols-length(hdr))...])
    Markdown.Table(name, tab(:l, cols))
end

function sumwidth(name::AbstractVector, cols::Integer)
    n = length(name)
    rows = (n + cols - 1) ÷ cols
    wi = 0
    for k = 1:cols
        r = (k-1)*rows+1:1:min(k*rows,n)
        if !isempty(r)
            wi += maximum(length.(name[r]))
        end
    end
    wi
end

vecvec(tab::Matrix) = [ tab[i,:] for i in 1:size(tab,1)]
tab(a, cols::Integer) = [a for i in 1:cols]

function buildnametable1(db::MatrixDatabase, hdr, p::Pattern, maxrow::Integer=0)
    dali = mdata.(Ref(db), mdlist(db, p))
    lt(a, b) = a.id < b.id
    sort!(dali, lt=lt)
    item(data::MatrixData) = string(data.id, " ", data.name)
    MD(buildnametable(hdr, item.(dali), maxrow))
end

MD(a...) = Markdown.MD([a...])

listnames(p::Pattern) = listnames(MATRIX_DB, p)
function listnames(db::MatrixDatabase, p::Pattern)
    li = mdlist(db, p)
    MD(buildnametable("list($(length(li)))", li))
end

##############################
# display database information
##############################

"""
    overview([db])

return formatted overview about matrices and groups in the collection.
"""
function overview(db::MatrixDatabase)
    # Print information strings
    LMARG = -10
    hdr_mat = Markdown.Header{3}("Currently loaded Matrices")

    dli(p) = mdata.(Ref(db), mdlist(db, p))
    ldall(p) = [length(mdlist(p & isloaded)), length(mdlist(p))]

    bmat = buildnametable1(db, "builtin(#)", :builtin, LMARG)
    umat = buildnametable1(db, "user(#)", :user, LMARG)
    spdir = buildnametable(["UFL  /  TAMU", "of"], ldall(sp(:)))
    mmdir = buildnametable(["MatrixMarket", "of"], ldall(mm(:)))

    hdr_groups = Markdown.Header{3}("Groups:")

    groups = buildnametable("Groups", listgroups(), LMARG)

    MD(([hdr_mat, bmat, umat, groups, spdir, mmdir]))
end

"""
    mdinfo()
Overview about matrices.
"""
mdinfo(db::MatrixDatabase=MATRIX_DB) = overview(db)

