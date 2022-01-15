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

@addgroup testgroup = ["rand1"]
@test mdlist(:testgroup) == ["rand1"]
@modifygroup testgroup = ["rand2"]
@test mdlist(:testgroup) == ["rand2"]
@rmgroup testgroup
@test_throws ArgumentError mdlist(:testgroup)

begin #Testing backward compatibility deprecation. Delete eventually.
    mydepot_warning = "MY_DEPOT_DIR custom code inclusion is deprecated: load custom generators by calling include_generator and reinitializing matrix depot at runtime. For more information, see: https://matrixdepotjl.readthedocs.io/en/latest/user.html. Duplicate warnings will be suppressed."

    mkpath(MatrixDepot.user_dir())
    open(joinpath(MatrixDepot.user_dir(), "group.jl"), "w") do f
        write(f, "usermatrixclass = Dict(\n);")
    end

    matrixgenerator = 
    """
    randorth(n) = Matrix(qr(randn(n,n)).Q)
    include_generator(FunctionName, "randorth", randorth)
    include_generator(Group, :random, randorth)
    """
    open(joinpath(MatrixDepot.user_dir(), "generator.jl"), "w") do f
        write(f, matrixgenerator)
    end

    @test_logs (:warn, mydepot_warning) min_level=Logging.Warn match_mode=:any MatrixDepot.init()
    n = rand(1:8)

    MatrixDepot.modgroup(:testgroup, ["rand1"])
    @test mdlist(:testgroup) == ["rand1"]

    @test matrixdepot("randorth", n) !== nothing
    @test mdinfo("randorth") !== nothing
    @test "randorth" in MatrixDepot.mdlist(:random)

    rm(joinpath(MatrixDepot.user_dir(), "group.jl"))
    rm(joinpath(MatrixDepot.user_dir(), "generator.jl"))
    @test_logs min_level=Logging.Warn MatrixDepot.init()
end

@test_throws ArgumentError include_generator(Group, :lkjasj, sin)
