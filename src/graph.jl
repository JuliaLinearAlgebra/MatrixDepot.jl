# adjacency matrices
#
# The code is adapted from CONTEST(Controlable TEST matrices) 
# by Alan Taylor and Professor Des Higham
# http://www.mathstat.strath.ac.uk/outreach/contest/

"""
Erdos-Renyi Random Graph
======================
Generate an adjacency matrix of an Edros-Renyi random graph:
an undirected graph is chosen uniformly at random from the set
of all symmetric graphs with `n` nodes and `m` edges.

*Input options:*

+ [type,] n, m: the dimension of the matrix is `n`. The number of
    1's in the matrix is `2*m`.

+ [type,] n: m = ceil(Int, n*log(n)/2).

*Groups:* ["sparse", "graph"]

*References:*

**P. Erdos and A. Renyi**, On Random Graphs, Publ. Math. Debrecen, 6, 1959, 
pp. 290-297
"""
function erdrey{T}(::Type{T}, n::Integer, m::Integer)
    nzeros = ceil(Int, 0.5*n*(n-1)*rand(m))
    v = zeros(Int, n)
    for count in 1:n
        v[count] = count*(count -one(Int))/2
    end
    
    is = zeros(Int, m)
    js = zeros(Int, m)
    
    for count in 1:m
        i = minimum(find(x -> x >=nzeros[count], v))
        j = nzeros[count] - (i-1)*(i-2)/2
        is[count] = i
        js[count] = j
    end
    A = sign(sparse([is;js], [js;is], ones(T, m*2), n, n))
    while nnz(A) != 2*m
        diff = round(Int, m - nnz(A)/2)
        is_new = zeros(Int, diff)
        js_new = zeros(Int, diff)
        for count = 1:diff
            idx = ceil(Int, 0.5*n*(n-1)*rand())
            is_new[count] = minimum(find(x -> x >= idx, v))
            js_new[count] = idx - (is_new[count] - 1)*(is_new[count] - 2)/2
        end
        append!(is, is_new)
        append!(js, js_new)
        s = ones(T, length(is))
        A = sign(sparse([is;js], [js;is], [s;s], n,n))
    end
    A
end
erdrey{T}(::Type{T}, n::Integer) = erdrey(T, n, ceil(Int, n*log(n)/2))
erdrey(arg...) = erdrey(Float64, arg...)

"""
Gilbert Random Graph
==================
Generate an adjecency matrix of a Gilbert random graph: an undirected graph
with pairs of nodes are connected with indepdent probability `p`.

*Input options:*

+ [type,] n, p: the dimension of the matrix is `n` and the probability that any two nodes
    are connected is `p`.

+ [type,] n: p = log(n)/n.

*Groups:* ["sparse", "graph"]

*References:*

**E.N. Gilbert**, Random Graphs, Ann. Math. Statist., 30, (1959) pp. 1141-1144.
"""
function gilbert{T}(::Type{T}, n::Integer, p::AbstractFloat)
    v = zeros(Int, n)
    for k = 1:n
        v[k] = round(Int, k*(k-1)/2)
    end
    
    is = zeros(Int, 0)
    js = zeros(Int, 0)
    
    w = zero(Int)

    w += one(Int) + floor(Int, log(1 - rand()) / log(1 - p))
    
    while w < n*(n-1)/2
        i = minimum(find(x -> x >= w, v))
        j = w - round(Int, (i -1)*(i - 2)/2)
        push!(is, i)
        push!(js, j)
        w += one(Int) + floor(Int, log(1 - rand()) / log(1-p))
    end
    
    s = ones(T, length(is))
    return sparse([is;js], [js;is], [s;s], n, n)
end
gilbert{T}(::Type{T}, n::Integer) = gilbert(T, n, log(n)/n)
gilbert(arg...) = gilbert(Float64, arg...)
