import MatrixDepot: include_generator, Group, FunctionName

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
MatrixDepot.init(ignoredb=true)

n = rand(1:8)
@test matrixdepot("randsym", n) !== nothing
@test mdinfo("randsym") != nothing
@test "randsym" in MatrixDepot.mdlist(:random)
@test "randsym" in MatrixDepot.mdlist(:symmetric)

@test_nowarn MatrixDepot.init()
n = rand(1:8)
@test matrixdepot("randsym", n) !== nothing
@test mdinfo("randsym") != nothing
@test "randsym" in MatrixDepot.mdlist(:random)
@test "randsym" in MatrixDepot.mdlist(:symmetric)

begin #Testing backward compatibility check deprecation. Delete eventually.
    mydepot_warning = "MY_DEPOT_DIR custom code inclusion is deprecated: load custom generators by calling include_generator and reinitializing matrix depot at runtime. For more information, see: https://matrixdepotjl.readthedocs.io/en/latest/user.html. Duplicate warnings will be suppressed."

    matrixdata = 
    """
    randorth(n) = Matrix(qr(randn(n,n)).Q)
    include_generator(FunctionName, "randorth", randorth)
    include_generator(Group, :random, randorth)
    """
    mkpath(MatrixDepot.user_dir())
    open(joinpath(MatrixDepot.user_dir(), "generator.jl"), "w") do f
        write(f, matrixdata)
    end

    @test_logs (:warn, mydepot_warning) MatrixDepot.init()
    n = rand(1:8)

    @test matrixdepot("randorth", n) !== nothing
    @test mdinfo("randorth") != nothing
    @test "randorth" in MatrixDepot.mdlist(:random)
end

@test_throws ArgumentError include_generator(Group, :lkjasj, sin)
