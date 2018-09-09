
using Markdown

"return information about all matrices selected by pattern"
mdinfo(p::Pattern) = mdinfo(MATRIX_DB, p)
function mdinfo(db::MatrixDatabase, p::Pattern)
    mdbuffer = Markdown.MD([])
    md = mdbuffer
    for name in list(p)
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

_mdheader(md::Markdown.MD, p, o) = isempty(md.content) ? (nothing, o) : _mdheader(md.content[1], md, o)
_mdheader(md::Markdown.Header, p, o) = (md, p)
_mdheader(md, p, o) = (nothing, o)

function mdinfo(data::GeneratedMatrixData)
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

function mdinfo(data::RemoteMatrixData)
    txt = mmreadcomment(matrixfile(data))
    txt = _repl(txt, r"^%-+$"m => "---", r"^%%" => "###### ", r"%+" => "* ")
    md = Markdown.parse(txt)
    insert!(md.content, 1, Markdown.Header{1}([data.name]))
    md
end

"""
format output list to multi- column table form adapted to display width.
"""
function buildnametable(name::AbstractVector, maxrow::Integer=0)
    maxrow = maxrow == 0 ? displaysize(stdout)[2] : maxrow
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
    insert!(name, 1, tab("", cols))    
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

function buildnametable1(db::MatrixDatabase, p::Pattern, maxrow::Integer=0)
    dali = mdata.(Ref(db), list(db, p))
    lt(a, b) = a.id < b.id
    sort!(dali, lt=lt)
    item(data::MatrixData) = string(data.id, " ", data.name)
    md(buildnametable(item.(dali), maxrow))
end

md(a...) = Markdown.MD([a...])

const mdlist = md ∘ buildnametable ∘ list

