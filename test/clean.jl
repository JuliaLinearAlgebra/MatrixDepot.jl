user_dir = joinpath(dirname(@__FILE__), "..", "myMatrixDepot")
dirdata = joinpath(dirname(@__FILE__),"..","data")
if isdir(user_dir)
    mv(user_dir, string(user_dir, ".saved"))
end
if isdir(dirdata)
    mv(user_dir, string(dirdata, ".saved"))
end

MatrixDepot.init()

