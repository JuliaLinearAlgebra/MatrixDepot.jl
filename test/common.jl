matrixdepot()
groups = ["symmetric", "inverse", "ill-cond", "pos-def", "eigen","sparse", "random", "regprob", "all", "data"]

for group in groups 
    matrixdepot(group)
end

@test_throws ArgumentError matrixdepot("something")

n = rand(1:55)
name = matrixdepot(n)
matrixdepot(name)

nlist = matrixdepot(1, 3, 4:20)
m = length(matrixdepot("all"))
@test_throws ArgumentError matrixdepot(m+1)

@addgroup newlist = matrixdepot(3:6, 20)

MatrixDepot.init()

println(matrixdepot("newlist"))

@rmgroup newlist

# testing the new API
#
@testset "load" begin
    # load all matrices loaded later in "test/download.jl"
    load([  "Bates/*"                
            "HB/1138_bus"                   
            "Harwell-Boeing/lanpro/nos5"    
            "Harwell-Boeing/smtape/bp___200"
            "Pajek/Journals"])
    @test length(list(:loaded)) == 7
end

@testset "list" begin

    LOC = length(list(:local))
@test LOC in [59, 60]   # depends on include_generator has worked 
    REM = length(list("*/*"))
@test REM in [2757, 2833]   # depends on whether ufl or tamu url has been used

@test list("") == []
@test list(1) == ["HB/1138_bus"]
@test list("#1") == ["HB/1138_bus"]
@test list("#M1") == ["Harwell-Boeing/psadmit/1138_bus"]
@test list("HB/1138_bus") == ["HB/1138_bus"]
@test length(list("HB/*")) == 292
@test list("Harwell-Boeing//") == ["Harwell-Boeing/*/* - (292)"]
@test list("*/hamm/*") == ["misc/hamm/add20", "misc/hamm/add32", "misc/hamm/memplus"]
@test list("*/hamm/*a*3?") == ["misc/hamm/add32"]
@test list("*/hamm/*a*3[123]") == ["misc/hamm/add32"]
@test list("*/hamm/*a*3[123]?") == []
@test length(list("*")) == LOC
@test length(list("*/*"))  == REM
@test length(list("*/*/*")) == 498
@test length(list("**")) == REM + 498 + LOC
@test list("*//*/*") == ["Harwell-Boeing/*/* - (292)", "NEP/*/* - (73)",
                           "SPARSKIT/*/* - (107)", "misc/*/* - (26)"]
@test list("//*") == ["/* - ($LOC)"]
@test list("//*/*") == ["/*/* - ($REM)"]
@test list("//*/*/*") == ["/*/*/* - (498)"]
@test list("HB/") == ["HB/* - (292)"]
@test length(list("Harwell-Boeing/*/*")) == 292
@test list(r".*ng/ma.*") == ["Harwell-Boeing/manteuffel/man_5976"]
@test list(2001:2002) == ["JGD_Groebner/c8_mat11_I", "JGD_Groebner/f855_mat9"]
@test length(list(2757:3000)) == REM - 2756
@test list(:xxx) == []
@test length(list(:remote)) == REM + 498
@test length(list(:local)) == LOC
@test length(list(:loaded)) + length(list(:unloaded)) == length(list(:remote))
@test length(list(:builtin)) + length(list(:user)) == length(list(:local))
@test length(list(:local)) + length(list(:remote)) == length(list(:all))
@test list(:all) == list("**")
@test length(list(:symmetric)) == 21
@test length(list(:illcond)) == 20
@test length(matrixdepot("ill-cond")) == 20

# intersections and unions
@test list((:posdef, :sparse)) == ["poisson", "wathen"]
@test length(list([:posdef,:sparse])) + length(list((:posdef,:sparse))) == length(list(:posdef)) + length(list(:sparse))


# predicates of remote matrices
@test length(list(issymmetric)) == 25 + LOC - 59

@test length(list(predm(n -> n < 10000))) == 6    # items with m < *
@test length(list(predn(n -> n < 10000))) == 5    # items with n < *
@test length(list(prednz(n -> n < 5000))) == 4    # items with nnz < *
@test length(list(predmn((m,n) -> m > n))) == 1   # items with m > n

end

