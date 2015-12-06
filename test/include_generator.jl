matrixdata = 
"""
"random symmetric matrix"
function randsym{T}(::Type{T}, n)
 A = zeros(T, n, n)
  for j = 1:n
      for i = j:n
          A[i,j] = randn()
      end
  end
  A = A + tril(A, -1)'
  return A
end
randsym(n) = randsym(Float64, n)
include_generator(FunctionName, "randsym", randsym)
include_generator(Group, "random", randsym)
include_generator(Group, "symmetric", randsym)
"""
user_dir = joinpath(dirname(@__FILE__), "..", "myMatrixDepot")

open(string(user_dir, "/generator.jl"), "w") do f
    write(f, matrixdata)
end
workspace()
using MatrixDepot
using Base.Test

n = rand(1:8)
A = matrixdepot("randsym", n)
matrixdepot("randsym")
@test "randsym" in matrixdepot("random")
@test "randsym" in matrixdepot("symmetric")

user_dir = joinpath(dirname(@__FILE__), "..", "myMatrixDepot")
rm(user_dir, recursive = true)
workspace()
using MatrixDepot
using Base.Test
