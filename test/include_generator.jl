import MatrixDepot: publish_user_generators, include_generator, Group, FunctionName
using Logging

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

# update the database
MatrixDepot.publish_user_generators()

n = rand(1:8)
@test matrixdepot("randsym", n) !== nothing
@test mdinfo("randsym") !== nothing
@test "randsym" in MatrixDepot.mdlist(:random)
@test "randsym" in MatrixDepot.mdlist(:symmetric)

@test_logs min_level=Logging.Warn MatrixDepot.init()
n = rand(1:8)
@test matrixdepot("randsym", n) !== nothing
@test mdinfo("randsym") !== nothing
@test mdinfo("randsym") == Base.Docs.doc(randsym)
@test "randsym" in MatrixDepot.mdlist(:random)
@test "randsym" in MatrixDepot.mdlist(:symmetric)

setgroup!(:testgroup, "randsy*")
@test mdlist(:testgroup) == ["randsym"]
setgroup!(:testgroup, ["*/1138_bus"])
@test mdlist(:testgroup) == ["HB/1138_bus"]
deletegroup!(:testgroup)
@test_throws ArgumentError mdlist(:testgroup)

@test_throws ArgumentError include_generator(Group, :lkjasj, sin)
