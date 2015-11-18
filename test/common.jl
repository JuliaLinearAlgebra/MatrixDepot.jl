matrixdepot()
groups = ["symmetric", "inverse", "ill-cond", "pos-def", "eigen","sparse", "random", "regprob", "all", "data"]

for group in groups 
    matrixdepot(group)
end

try 
    matrixdepot("something")
catch ArgumentError
    println("matrixdepot: no information for `something`")
end

n = rand(1:55)
name = matrixdepot(n)
matrixdepot(name)

list = matrixdepot(1, 3, 4:20)

@addgroup newlist = matrixdepot(3:6, 20)

workspace()
using MatrixDepot
using Base.Test

@rmgroup newlist

