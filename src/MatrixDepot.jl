module MatrixDepot
export matrixdepot

include("higham.jl") #Higham Test matrices

function matrixdepot{T}(name::String, ::Type{T}, m::Int, n::Int)
    # name is matrix name
    # m is the number of rows
    # n is the number of columns
    return matrixdict[name](T, m, n)
end
matrixdepot(name::String, m::Int, n::Int) = matrixdepot(name, Float64, m, n)

function matrixdepot{T}(name::String, ::Type{T}, n::Int)
    # name is the matrix name
    # n is the dimension of the matrix (square)
    return matrixdict[name](T, n)
end
matrixdepot(name::String, n::Int) = matrixdepot(name, Float64, n)

function matrixdepot{T}(name::String, x::Vector{T}, y::Vector{T})
    # name: matrix name
    # x,y : vectors
    return matrixdict[name](x,y)
end
matrixdepot{T}(name::String, x::Vector{T}) = matrixdepot(name, x, x)

function matrixdepot{T}(name::String, x::Vector{T}, n::Int)
    # name: matrix name
    # x: a vector
    # the column dimension of the matrix
    return matrixdic[name](x,n)
end

function matrixdepot(name::String)
    # name is the matrix name
    if name in keys(matrixinfo)
        println(matrixinfo[name])
    elseif name in keys(matrixclass)
        return matrixclass[name]
    else
        println("Your matrix or class is not included in Matrix Depot.")
    end
end

end # end module
