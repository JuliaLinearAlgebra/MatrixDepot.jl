# Test matrices for regularization methods from Hansen's
# Regularization toolbox

immutable RegProb{T}
    A::Matrix{T}  # matrix of interest
    b::Vector{T}  # right-hand side
    x::Vector{T}  # the solution to Ax = b
end

# compute a oscillating matrix from a given signular value spectrum
# using the idea proposed by Per Christian Hansen, "Test matrices for
# regularization methods. 
function oscillate{T}(Σ::Vector{T})
    n = length(Σ)
    dv = rand(T, 1, n)[:] + eps(T)
    ev = rand(T, 1, n-1)[:] + eps(T)
    B = Bidiagonal(dv, ev, true)
    U, S, V = svd(B)
    return U*diagm(Σ)*U'
end
        
#  mode = 1: one large singular value.
#  mode = 2: one small singular value.
#  mode = 3: geometrically distributed singular values.
#  mode = 4: arithmetrically distributed singular values.
#  mode = 5: random singular values with  unif. dist. logarithm.
#  κ = sqrt(1/eps(T)) is the condition number of the matrix.
function oscillate{T}(::Type{T}, n::Int, mode::Int)
    κ = sqrt(1/eps(T))
    if mode == 3
        factor = κ^(-1/(n-1))
        Σ = factor.^[0:n-1;]
    elseif mode == 4
        Σ = ones(T, n) - T[0:n-1;]/(n-1)*(1 - 1/κ)
    elseif mode == 2
        Σ = ones(T, n)
        Σ[n] = one(T)/κ
    elseif mode == 1
        Σ = ones(n)./κ
        Σ[1] = one(T)
    else
        error("invalid mode value.")
    end
    return oscillate(Σ)
end
oscillate{T}(::Type{T}, n::Int) = oscillate(T, n, 3)

