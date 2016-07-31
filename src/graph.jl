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
function gilbert{T}(::Type{T}, n::Integer)
    if n == 1
        return gilbert(T, n, 0.2) 
    else
        return gilbert(T, n, log(n)/n)
    end
end
gilbert(arg...) = gilbert(Float64, arg...)


# utility function
# shortcuts: randomly add entries (shortcuts) to a matrix.
function shortcuts{T}(A::SparseMatrixCSC{T}, p::Real)
    n, = size(A)
    Ihat = find(x -> x <= p, rand(n))
    Jhat = ceil(Int, n*rand(size(Ihat)))
    Ehat = ones(T, size(Ihat))
    
    # an edge
    for (i, ele) in enumerate(Ihat)
        if ele == Jhat[i]
            Ehat[i] = zero(T)
        end
    end

    is, js, es = findnz(A)
    
    return sparse([is; Ihat; Jhat], [js; Jhat; Ihat], [es; Ehat; Ehat], n, n)
end

"""
Small World Network
=================
Generate an adjacency matrix for a small world network. We model
it by adding shortcuts to a kth nearest neighbour ring network 
(nodes i and j are connected iff |i -j| <= k or |n - |i -j|| <= k.) with n nodes. 

*Input options:* 

+ [type,] n, k, p: the dimension of the matrix is `n`. The number of 
    nearest-neighbours to connect is `k`. The probability of adding a shortcut
    in a given row is `p`.

+ [type,] n: `k = 2` and `p = 0.1`.

*References:*

**D. J. Watts and S. H. Strogatz**, Collective Dynamics of Small World Networks,
Nature 393 (1998), pp. 440-442.
"""
function smallworld{T}(::Type{T}, n::Integer, k::Integer, p::Real)
    twok = 2*k
    is = zeros(Int, 2*k*n)
    js = zeros(Int, 2*k*n)
    es = zeros(T, 2*k*n)

    for count = 1:n
        is[(count-1)*twok+1 : count*twok] = count*ones(Int, twok)
        js[(count-1)*twok+1 : count*twok] = mod([count : count+k-1; 
                                                                               n-k+count-1 : n+count-2], n) + 1
        es[(count-1)*twok+1 : count*twok] = ones(T, twok)
    end

    A = sparse(is, js, es, n, n)
    return sign(shortcuts(A, p))
end
smallworld{T}(::Type{T}, n::Integer) = smallworld(T, n, 2, 0.1)
smallworld(arg...) = smallworld(Float64, arg...)
