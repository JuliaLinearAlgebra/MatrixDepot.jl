dirdata = joinpath(dirname(@__FILE__),"..","data")
if isdir(string(dirdata, '/', "uf"))
    rm(string(dirdata, '/', "uf"), recursive = true)
end
if isdir(string(dirdata, '/', "mm"))
    rm(string(dirdata, '/', "mm"), recursive = true)
end

#download
# MatrixDepot.update() # already done when initializing directory
# uf
matrixdepot("1138_bus", :get)
matrixdepot("HB/1138_bus", :get)
matrixdepot("Pajek/Journals", :get)
# Matrix Markt
matrixdepot("Harwell-Boeing/smtape/bp___200", :get)
 
# read data
A = matrixdepot("HB/1138_bus", :r)
matrixdepot("HB/1138_bus")
matrixdepot("1138_bus", :s) #search
#B = matrixdepot("Harwell-Boeing/psadmit/662_bus", :read)
#matrixdepot("Harwell-Boeing/psadmit/662_bus")
r = matrixdepot("Pajek/Journals", :read, meta = true)
display(r["Journals"])
# MatrixDepot.update() # already done

matrixdepot("Bates/*", :get)
B = matrixdepot("Bates/Chem97Zt", :r)
matrixdepot("Bates/Chem97Zt")
matrixdepot()
@test_throws ArgumentError matrixdepot("HB/662_bus", :k)

matrixdepot("epb0", :get)

# matrix market
matrixdepot("Harwell-Boeing/lanpro/nos5", :get)
matrixdepot()

# an example with rhs and solution
data = mdopen("DRIVCAV/cavity01")
@test MatrixDepot.loadmatrix(data) === nothing
@test size(matrix(data)) == (317, 317)
@test size(rhs(data), 1) == 317
@test size(solution(data), 1) == 317

# an example loading a txt file
data = mdopen("Pajek/Journals")
@test length(metareader(data, "Journals_nodename.txt")) > 100


# rm data

#=
rm(string(dirdata, '/', "uf_matrices.html"))
rm(string(dirdata, '/', "mm_matrices.html"))
rm(string(dirdata, '/', "uf"), recursive = true)
rm(string(dirdata, '/', "mm"), recursive = true)
=#

