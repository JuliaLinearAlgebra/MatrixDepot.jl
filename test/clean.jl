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
        ufm = "uf_matrices.html"
        cp(abspath(target, ufm), abspath(path, ufm))
        mmm = "mm_matrices.html"
        cp(abspath(target, mmm), abspath(path, mmm))
    end
    target
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
