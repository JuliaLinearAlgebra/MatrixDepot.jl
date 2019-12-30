mdinfo()
groups = ["symmetric", "inverse", "illcond", "posdef", "eigen","sparse", "random", "regprob", "all"]

for group in groups 
    listnames(Symbol(group))
end

@test_throws DataError matrixdepot("something")

n = rand(1:55)
name = mdlist(builtin(n))
@test mdinfo(name) !== nothing

nlist = mdlist(builtin([1, 3]) | builtin(4:20))
m = length(mdlist("**"))
@test isempty(mdinfo(builtin(m+1)))

@addgroup newlist = mdlist(builtin(3:6) | builtin(20))

MatrixDepot.init(ignoredb=true)

@test mdlist(:newlist) == mdlist(builtin(3:6,20))

@rmgroup newlist

# testing the new API
#
# it is assumed that all matrices generated during "test/dowlnload.jl"
# and the ibuiltin- and user-defined "randsym" but no others are available
#

@testset "mdlist" begin

REM = length(mdlist("*/*"))
@test REM in [2757, 2856]   # depends on whether UFL or TAMU url has been used
@test length(mdlist(:builtin)) == 59
@test length(mdlist(:user)) in [0, 1]

@test mdlist("") == []
@test mdlist("HB/1138_bus") == ["HB/1138_bus"]
@test mdlist(sp(1)) == ["HB/1138_bus"]
@test mdlist(mm(1)) == ["Harwell-Boeing/psadmit/1138_bus"]
@test sort(mdlist(sp(1:3000))) == mdlist("*/*")
@test sort(mdlist(mm(1:3000))) == mdlist("*/*/*")
@test mdlist(builtin(:)) == mdlist(isbuiltin)
@test mdlist(user(:)) == mdlist(isuser)
@test mdlist("*") == mdlist(islocal)
@test length(mdlist("HB/*")) == 292
@test mdlist("HB**") == mdlist("HB/**")
@test_throws ArgumentError  listdir("*/*")
@test listdir("Harwell-Boeing//") == ["Harwell-Boeing/*/* - (292)"]
@test mdlist("*/hamm/*") == ["misc/hamm/add20", "misc/hamm/add32", "misc/hamm/memplus"]
@test mdlist("*/hamm/*a*3?") == ["misc/hamm/add32"]
@test mdlist("*/hamm/*a*3[123]") == ["misc/hamm/add32"]
@test mdlist("*/hamm/*a*3[123]?") == []
@test length(mdlist("*/*/*")) == 498

@test listdir("*//*/*") == ["Harwell-Boeing/*/* - (292)", "NEP/*/* - (73)",
                           "SPARSKIT/*/* - (107)", "misc/*/* - (26)"]
@test listdir("//*") == ["/* - ($(length(mdlist(:local))))"]
@test listdir("//*/*") == ["/*/* - ($REM)"]
@test listdir("//*/*/*") == ["/*/*/* - (498)"]
@test listdir("HB/") == ["HB/* - (292)"]
@test length(mdlist("Harwell-Boeing/*/*")) == 292
@test mdlist(r".*ng/ma.*") == ["Harwell-Boeing/manteuffel/man_5976"]
@test mdlist(sp(2001:2002)) == ["JGD_Groebner/c8_mat11_I", "JGD_Groebner/f855_mat9"]
@test length(mdlist(sp(2757:3000))) == REM - 2756
@test_throws ArgumentError mdlist(:xxx)
@test length(mdlist(isremote)) == REM + 498
@test length(mdlist(isloaded)) + length(mdlist(isunloaded)) == length(mdlist(isremote))
@test length(mdlist(isbuiltin)) + length(mdlist(isuser)) == length(mdlist(islocal))
@test length(mdlist(islocal)) + length(mdlist(isremote)) == length(mdlist(:all))
@test length(mdlist("**")) == length(mdlist("*")) + length(mdlist("*/*")) + length(mdlist("*/*/*"))
@test mdlist(:all) == mdlist("**")
@test length(mdlist(:symmetric)) == 22
@test length(mdlist(:illcond)) == 20

# intersections and unions
@test mdlist((:posdef, :sparse)) == ["poisson", "wathen"]
@test length(mdlist([:posdef,:sparse])) + length(mdlist((:posdef,:sparse))) == length(mdlist(:posdef)) + length(mdlist(:sparse))


# predicates of remote and local matrices
@test length(mdlist(issymmetric)) == 30

@test length(mdlist(@pred(m < 10000))) == 1672   # items with m < *
@test length(mdlist(@pred(m < 10000) & isloaded)) == 10   # items with m < *
@test length(mdlist(@pred(n < 10000))) == 1612   # items with n < *
@test length(mdlist(@pred(n < 10000) &  isloaded)) == 9    # items with n < *
@test length(mdlist(@pred(nnz < 5000))) == 600    # items with nnz < *
@test length(mdlist(@pred(nnz < 5000) & isloaded)) == 2    # items with nnz < *
@test length(mdlist(@pred(m > n))) == 178   # items with m > n
@test length(mdlist(@pred(m > n) & isloaded)) == 0   # items with m > n
@test length(mdlist(@pred(occursin("Power", kind)))) == 70
@test length(mdlist(@pred(0 < date <= 1971))) == 5
@test length(mdlist(@pred(date >= 2016) & @pred(0< date <=2016))) == 2
@test length(mdlist(@pred(date == 0) & sp(:))) == 42
@test length(mdlist(:symmetric & "kahan")) == 0
@test length(mdlist(:symmetric & "hankel")) == 1

@test mdlist(:local) == mdlist(islocal)
@test mdlist(:builtin) == mdlist(isbuiltin)
@test mdlist(:user) == mdlist(isuser)

@test mdlist(:symmetric) == mdlist(issymmetric & islocal)
end

@testset "logical" begin
# for the boolean syntax
@test mdlist(~islocal) == mdlist(isremote)
@test length(mdlist(isloaded & issymmetric)) == 7
@test length(mdlist(isloaded & ~issymmetric)) == 4
@test length(mdlist(isloaded & ~issymmetric | isuser & issymmetric)) == 5
@test length(mdlist(!islocal & issymmetric | isuser & issymmetric)) == 9
@test mdlist(islocal & ~isbuiltin) == mdlist(isuser)
@test mdlist(islocal & ~isbuiltin) == mdlist(isuser)
@test "a" & "b" === ("a", "b")
@test ~"a" & "b" === (~"a", "b")
@test ~"a" * "b" === ~"ab"
@test ~"ab" === ~'a' * 'b'
@test ~~islocal === islocal
@test mdlist(~"*a*" & ~ "*e*") == mdlist(~["*a*", "*e*"])
@test mdlist(()) == mdlist(:all)
@test "a" & r"b" == ("a", r"b")
@test "a" & r"b" & "c" == ("a", r"b", "c")
@test "a" | r"b" == ["a", r"b"]
@test "a" | r"b" | "c" == ["a", r"b", "c"]
@test mdlist(islocal & ~("*a*" | "*e*")) == mdlist(:local & ~"*a*" & ~ "*e*")
@test "a" & ( "b" & "c" ) == ("a", "b", "c")
@test "a" | ( "b" | "c" ) == ["a", "b", "c"]
@test_throws ArgumentError mdlist([] & :invalid_group_name)

@test MatrixDepot.mdlist(builtin(10, 1:7, 3)) == ["baart", "binomial", "blur", "cauchy",
                                                "chebspec", "chow", "circul", "deriv2"] 

@test MatrixDepot.fname(sin) == "unknown-function"
@test_throws MethodError mdopen("baart", 10, 11)
end

@testset "properties" begin

md = mdopen("gravity", 10)
io = IOBuffer()
@test show(io, md) === nothing
@test propertynames(md) == [:A, :m, :n, :nnz, :dnz]
@test md.A isa AbstractMatrix
@test md.nnz == count(md.A .!= 0)
@test md.dnz == 0
@test_throws DataError md.x

md = mdopen("gravity", 10, false)
io = IOBuffer()
@test show(io, md) === nothing
@test propertynames(md) == [:A, :b, :x, :m, :n, :nnz, :dnz]
@test md.A isa AbstractMatrix
@test md.b isa AbstractVector
@test md.x isa AbstractVector

md = mdopen("HB/well1033")
io = IOBuffer()
@test show(io, md) === nothing
@test propertynames(md) == [:A, :b, :m, :n, :nnz, :dnz]
@test md.A isa AbstractMatrix
@test md.b isa AbstractMatrix
@test_throws DataError md.x
@test_throws DataError md.invalidname
@test_throws ErrorException md.data.invalidname
@test md.dnz == md.data.header.nnz

@test length(mdlist(keyword("graph" & ("butterfly" | "network")) & isloaded)) == 2
@test length(mdlist(hasdata(:A & (:b | :x)))) == 2

end

