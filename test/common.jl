matrixdepot()
groups = ["symmetric", "inverse", "ill-cond", "pos-def", "eigen","sparse", "random", "regprob", "all", "data"]

for group in groups 
    matrixdepot(group)
end

@test_throws ArgumentError matrixdepot("something")

n = rand(1:55)
name = matrixdepot(n)
matrixdepot(name)

nlist = matrixdepot(1, 3, 4:20)
m = length(matrixdepot("all"))
@test_throws ArgumentError matrixdepot(m+1)

@addgroup newlist = matrixdepot(3:6, 20)

MatrixDepot.init()

println(matrixdepot("newlist"))


@rmgroup newlist

