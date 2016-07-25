# adjacency matrices for graphs


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
        diff = m - nnz(A)/2
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
