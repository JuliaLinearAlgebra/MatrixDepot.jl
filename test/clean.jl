# setup clean and empty directories

user_dir = abspath(dirname(@__FILE__), "..", "myMatrixDepot")
dirdata = abspath(dirname(@__FILE__),"..","data")

function savetarget(path::AbstractString)
    base = basename(path)
    dir = dirname(path)
    target = abspath(dir, string(base, ".saved"))
    while isdir(target)
        target = abspath(target, string(base, ".saved"))
    end
    if isdir(path)
        mv(path, target)
    end
end

savetarget(user_dir)
savetarget(dirdata)

# that will lownload the index files and initialize internal data
MatrixDepot.init()

