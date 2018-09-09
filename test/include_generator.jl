matrixdata = 
"""
"random symmetric matrix"
function randsym(::Type{T}, n) where T
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
include_generator(Group, :random, randsym)
include_generator(Group, :symmetric, randsym)
"""
user_dir = abspath(dirname(@__FILE__), "..", "myMatrixDepot")

open(joinpath(user_dir, "generator.jl"), "w") do f
    write(f, matrixdata)
end

# load the just written user file
MatrixDepot.init()

n = rand(1:8)
@test matrixdepot("randsym", n) !== nothing
@test mdinfo("randsym") != nothing
@test "randsym" in MatrixDepot.list(:random)
@test "randsym" in MatrixDepot.list(:symmetric)

import MatrixDepot: include_generator, Group
@test_throws ArgumentError include_generator(Group, :lkjasj, sin)

