# setup clean and empty directories

user_dir = abspath(dirname(@__FILE__), "..", "myMatrixDepot")
data_dir = abspath(dirname(@__FILE__),"..","data")

function save_target(path::AbstractString)
    base = basename(path)
    dir = dirname(path)
    target = abspath(dir, string(base, ".saved"))
    while isdir(target)
        target = abspath(target, string(base, ".saved"))
    end
    if isdir(path)
        mv(path, target)
    end
    if basename(path) == "data"
        mkpath(path)
        cp_if("uf_matrices.html", target, path)
        cp_if("mm_matrices.html", target, path)
        cp_if("ss_index.mat", target, path)
    end
    target
end

function cp_if(name::AbstractString, from::AbstractString, to::AbstractString)
    fromname = abspath(from, name)
    toname = abspath(to, name)
    if isfile(fromname) && !MatrixDepot.URL_REDIRECT[]
        cp(fromname, toname)
    end
end

function revert_target(saved::AbstractString, orig::AbstractString)
    savetest = abspath(dirname(orig), string(basename(orig), ".test"))
    if isdir(orig)
        if isdir(savetest)
            mv(savetest, joinpath(orig, basename(savetest)))
        end
        mv(orig, savetest)
    end
    if isdir(saved)
        mv(saved, orig)
    end
end
