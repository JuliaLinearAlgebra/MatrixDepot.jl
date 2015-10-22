dirdata = joinpath(Pkg.dir("MatrixDepot"), "data")
if isdir(string(dirdata, '/', "uf"))
    rm(string(dirdata, '/', "uf"), recursive = true)
end
if isdir(string(dirdata, '/', "mm"))
    rm(string(dirdata, '/', "mm"), recursive = true)
end

#download
MatrixDepot.update()
matrixdepot("HB/1138_bus", :get)
#matrixdepot("Harwell-Boeing/psadmit/662_bus", :get)
# read data
A = matrixdepot("HB/1138_bus", :r)
matrixdepot("HB/1138_bus")
#B = matrixdepot("Harwell-Boeing/psadmit/662_bus", :read)
#matrixdepot("Harwell-Boeing/psadmit/662_bus")

# rm data

rm(string(dirdata, '/', "uf_matrices.html"))
rm(string(dirdata, '/', "mm_matrices.html"))
rm(string(dirdata, '/', "uf"), recursive = true)
rm(string(dirdata, '/', "mm"), recursive = true)
