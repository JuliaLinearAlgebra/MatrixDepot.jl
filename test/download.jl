#download
MatrixDepot.update()
MatrixDepot.get("HB/1138_bus")

# read data
A = matrixdepot("HB/1138_bus", :r)
matrixdepot("HB/1138_bus")

# rm data
dirdata = joinpath(Pkg.dir("MatrixDepot"), "data")
rm(string(dirdata, '/', "uf_matrices.html"))
rm(string(dirdata, '/', "uf"), recursive = true)
