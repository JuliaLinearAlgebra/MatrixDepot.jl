# Test matrices for regularization methods from Hansen's
# Regularization toolbox
#
# Per Christian Hansen, Test matrices for regularization methods
# SIAM J. SCI. COMPUT Vol 16, No2, pp 506-512, (1995)


immutable RegProb{T}
    A::Matrix{T}  # matrix of interest
    b::Vector{T}  # right-hand side
    x::Vector{T}  # the solution to Ax = b
end


function baart{T}(::Type{T}, n::Int)

end
