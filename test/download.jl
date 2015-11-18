dirdata = joinpath(dirname(@__FILE__),"..","data")
if isdir(string(dirdata, '/', "uf"))
    rm(string(dirdata, '/', "uf"), recursive = true)
end
if isdir(string(dirdata, '/', "mm"))
    rm(string(dirdata, '/', "mm"), recursive = true)
end

#download
MatrixDepot.update()
# uf
matrixdepot("HB/1138_bus", :get)
# Matrix Markt
#matrixdepot("Harwell-Boeing/smtape/bp___200", :get)
 
# read data
A = matrixdepot("HB/1138_bus", :r)
matrixdepot("HB/1138_bus")
matrixdepot("1138_bus", :s) #search
#B = matrixdepot("Harwell-Boeing/psadmit/662_bus", :read)
#matrixdepot("Harwell-Boeing/psadmit/662_bus")
MatrixDepot.update()

matrixdepot("Bates/*", :get)
B = matrixdepot("Bates/Chem97Zt", :r)
matrixdepot("Bates/Chem97Zt")

try
    matrixdepot("HB/662_bus", :g)
catch ArgumentError
    println("download error: unknown argument :g")
end
# rm data


rm(string(dirdata, '/', "uf_matrices.html"))
rm(string(dirdata, '/', "mm_matrices.html"))
rm(string(dirdata, '/', "uf"), recursive = true)
rm(string(dirdata, '/', "mm"), recursive = true)

