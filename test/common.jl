matrixdepot()
groups = ["symmetric", "inverse", "ill-cond", "pos-def", "eigen","sparse", "random", "regprob", "all", "data"]

for group in groups 
    matrixdepot(group)
end

@test_throws DataError matrixdepot("something")

n = rand(1:55)
name = matrixdepot(n)
matrixdepot(name)

nlist = matrixdepot(1, 3, 4:20)
m = length(matrixdepot("all"))
@test_throws DataError matrixdepot(m+1)

@addgroup newlist = matrixdepot(3:6, 20)

MatrixDepot.init()

println(matrixdepot("newlist"))

@rmgroup newlist

# testing the new API
#
# it is assumed that all matrices generated during "test/dowlnload.jl"
# and the ibuiltin- and user-defined "randsym" but no others are available
#

@testset "list" begin

REM = length(list("*/*"))
@test REM in [2757, 2833]   # depends on whether ufl or tamu url has been used
@test length(list(:builtin)) == 59
@test length(list(:user)) in [0, 1]

@test list("") == []
@test list("HB/1138_bus") == ["HB/1138_bus"]
@test list(1) == ["HB/1138_bus"]
@test list("#1") == ["HB/1138_bus"]
@test list("#M1") == ["Harwell-Boeing/psadmit/1138_bus"]
@test list("#[123456789]*") == list("*/*")
@test list("#M*") == list("*/*/*")
@test list("#B*") == list(:builtin)
@test list("#U*") == list(:user)
@test list("*") == list(:local)
@test length(list("HB/*")) == 292
@test list("Harwell-Boeing//") == ["Harwell-Boeing/*/* - (292)"]
@test list("*/hamm/*") == ["misc/hamm/add20", "misc/hamm/add32", "misc/hamm/memplus"]
@test list("*/hamm/*a*3?") == ["misc/hamm/add32"]
@test list("*/hamm/*a*3[123]") == ["misc/hamm/add32"]
@test list("*/hamm/*a*3[123]?") == []
@test length(list("*/*/*")) == 498
@test list("*//*/*") == ["Harwell-Boeing/*/* - (292)", "NEP/*/* - (73)",
                           "SPARSKIT/*/* - (107)", "misc/*/* - (26)"]
@test list("//*") == ["/* - ($(length(list(:local))))"]
@test list("//*/*") == ["/*/* - ($REM)"]
@test list("//*/*/*") == ["/*/*/* - (498)"]
@test list("HB/") == ["HB/* - (292)"]
@test length(list("Harwell-Boeing/*/*")) == 292
@test list(r".*ng/ma.*") == ["Harwell-Boeing/manteuffel/man_5976"]
@test list(2001:2002) == ["JGD_Groebner/c8_mat11_I", "JGD_Groebner/f855_mat9"]
@test length(list(2757:3000)) == REM - 2756
@test list(:xxx) == []
@test length(list(:remote)) == REM + 498
@test length(list(:loaded)) + length(list(:unloaded)) == length(list(:remote))
@test length(list(:builtin)) + length(list(:user)) == length(list(:local))
@test length(list(:local)) + length(list(:remote)) == length(list(:all))
@test length(list("**")) == length(list("*")) + length(list("*/*")) + length(list("*/*/*"))
@test list(:all) == list("**")
@test length(list(:symmetric)) == 23
@test length(list(:illcond)) == 20
@test length(matrixdepot("ill-cond")) == 20

# intersections and unions
@test list((:posdef, :sparse)) == ["poisson", "wathen"]
@test length(list([:posdef,:sparse])) + length(list((:posdef,:sparse))) == length(list(:posdef)) + length(list(:sparse))


# predicates of remote and local matrices
@test length(list(issymmetric)) == 27

@test length(list(predm(n -> n < 10000))) == 9    # items with m < *
@test length(list(predn(n -> n < 10000))) == 8    # items with n < *
@test length(list(prednz(n -> n < 5000))) == 5    # items with nnz < *
@test length(list(predmn((m,n) -> m > n))) == 0   # items with m > n

# for the boolean syntax
@test list(!islocal) == list(isremote)
@test length(list(isloaded & issymmetric)) == 5
@test length(list(isloaded & !issymmetric)) == 4
@test length(list(isloaded & !issymmetric | isuser & issymmetric)) == 5
@test length(list(!islocal & issymmetric | isuser & issymmetric)) == 6
@test list(islocal & !isbuiltin) == list(isuser)

end

