user_dir = joinpath(dirname(@__FILE__), "..", "myMatrixDepot")
dirdata = joinpath(dirname(@__FILE__),"..","data")
 if isdir(user_dir)
     rm(user_dir, recursive = true)
end
if isdir(string(dirdata, '/', "uf"))
    rm(string(dirdata, '/', "uf"), recursive = true)
end
if isdir(string(dirdata, '/', "mm"))
    rm(string(dirdata, '/', "mm"), recursive = true)
end
workspace()
using MatrixDepot
using Base.Test
