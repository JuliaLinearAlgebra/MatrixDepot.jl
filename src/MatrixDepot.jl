module MatrixDepot
export matrixdepot

include("higham.jl") #Higham Test matrices

function matrixdepot()
    # Print information strings 
    println("Matrices in the Collection:")
    println("---------------------------")
    println()
    for prop in keys(matrixclass)
        println("Property: ", prop)
        println("-------------------")
        count = 0
        for mat in matrixclass[prop]
        print(mat, " matrix, ")
            count += 1
            if count > 3
                count = 0
                print("\n")
            end
        end
        println()
        println()
    end
end

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

function matrixdepot(name::String, n::Int)
    # name is the matrix name
    # n is the dimesion of the matrix
    # magic square is an exception: Int Array by default.
    if name == "magic"
        matrixdepot(name, Int, n)
    else
        matrixdepot(name, Float64, n)
    end
end

function matrixdepot{T}(name::String, ::Type{T}, n::Int, alpha, beta)
    # name is the matrix name
    # n is the dimension of the matrix
    # alpha, beta are scalars
    return matrixdict[name](T, n, alpha, beta)
end
matrixdepot(name::String, n::Int, alpha, beta) = matrixdepot(name, Float64, n, alpha, beta)

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
