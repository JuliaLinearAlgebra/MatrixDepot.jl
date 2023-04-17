#download

import MatrixDepot: load, loadinfo, loadsvd

# for test coverage
MatrixDepot.update()

# verify that no data are downloaded (empty directory)
@test mdlist(isloaded) == []

# Suite Sparse
@test loadinfo("*/1138_bus") == 1 # load only header
@test load("**/1138_bus") == 2 # loaded both versions (sp and mm)
@test load("HB/1138_bus") == 0 # count actually loaded
@test load("Pajek/Journals") == 1
@test load("wing") == 0 # do not try to load local matrix

# Matrix Market
@test load("Harwell-Boeing/smtape/bp___200") == 1

@test mdlist(isloaded) == ["HB/1138_bus", "Harwell-Boeing/psadmit/1138_bus",
                        "Harwell-Boeing/smtape/bp___200", "Pajek/Journals"]

# read data
@test matrixdepot("HB/1138_bus") !== nothing
@test mdlist("**/1138_bus") == ["HB/1138_bus", "Harwell-Boeing/psadmit/1138_bus"]
@test string(mdinfo("Harwell-Boeing/psadmit/662_bus")) ==
        "# Harwell-Boeing/psadmit/662_bus\n\n" *
        "###### MatrixMarket matrix coordinate real symmetric\n\n662 662 1568\n"

@test loadinfo(MatrixDepot.mdata("HB/1138_bus")) == 0
@test loadinfo(MatrixDepot.mdata("Bates/Chem97Zt")) == 1
@test length(string(mdinfo("Bates/Chem97Zt"))) == 447

# metadata access with old interface
@test MatrixDepot.metareader(mdopen("Pajek/Journals"), "Journals.mtx") !== nothing
@test_throws DataError MatrixDepot.metareader(mdopen("*/1138_bus"), "1138_bus_b")
@test_throws DataError MatrixDepot.metareader(mdopen("that is nothing"))

@test load("Bates/C*") == 2
@test mdinfo("Bates/Chem97Zt") !== nothing
@test matrixdepot("Bates/Chem97Zt") !== nothing
@test "Bates/Chem97Zt" in mdlist(isloaded)
@test length(mdlist(isloaded)) == 6

@test load("**/epb0") == 1

@test loadsvd("*/1138_bus") == 1

# matrix market
@test load("Harwell-Boeing/lanpro/nos5") == 1
@test length(mdlist(isloaded)) == 8
@test mdopen("Bai/dwg961b") !== nothing
mdesc = mdopen("Bai/dwg961b")
@test iscomplex(mdesc.data)
@test issymmetric(mdesc.data)
@test !ishermitian(mdesc.data)
@test !isboolean(mdesc.data)
@test metasymbols(mdesc) == [:A]
@test mdesc.A == matrixdepot("Bai/dwg961b")
mdesc = mdopen("Harwell-Boeing/cegb/cegb2802")
@test isloaded(mdesc.data) == false
@test_throws DataError mdesc.A

# an example with rhs and solution
mdesc = mdopen("DRIVCAV/cavity14")
A1 = mdesc.A
@test size(A1) == (2597, 2597)
A2 = mdesc.A
@test A1 === A2 # cache should be used
@test size(mdesc.b, 1) == 2597
@test size(mdesc.x, 1) == 2597
@test_throws DataError mdesc["invlid"]
@test_throws DataError mdesc.y

# read a format array file
@test MatrixDepot.mmreadheader(abspath(MatrixDepot.matrixfile(mdesc.data),
                                       "..", string("cavity14_b.mtx"))) !== nothing

mdesc = mdopen("blur", 10, false)
@test mdesc.A isa AbstractMatrix
@test mdesc.b isa AbstractVector
@test mdesc.x isa AbstractVector
@test_throws DataError mdesc["invlid"]
@test_throws DataError mdesc.y

@test_throws DataError mdopen(sp(1)).b # no rhs data for this example

mdesc = mdopen("*/bfly")
@test mdesc["Gname_01.txt"] == "BFLY3\n"
@test mdesc("Gname_01.txt") == "BFLY3\n"
@test size(mdesc.G_06) == (2048, 2048)
@test mdesc.G_06 === mdesc["G_06"]
@test mdesc.G_06 === mdesc[Symbol("G_06.mtx")]
fn = joinpath(dirname(MatrixDepot.matrixfile(mdesc.data)), "bfly_Gname_01.txt")
@test_throws DataError MatrixDepot.mmreadheader(fn)

# read from a pipeline
io = IOBuffer("%%matrixmarket matrix array real general\n1 1\n2.5\n")
@test MatrixDepot.mmread(io) == reshape([2.5], 1, 1)

# read from temporary file
mktemp() do path, io
    println(io, "%%MatrixMarket Matrix arRAy real GeneraL")
    println(io, "% author: willi")
    println(io)
    println(io, "% notes:")
    println(io, "%second Line")
    println(io); println(io, "    ")
    println(io, "24 2")
    close(io)
    @test MatrixDepot.mmreadcomment(path) ==
    "%%MatrixMarket Matrix arRAy real GeneraL\n% author: willi\n\n" *
    "% notes:\n%second Line\n\n    \n24 2\n"
    @test MatrixDepot.mmreadheader(path) ==
    Dict{Symbol,Any}(:symmetry=>"general",:m=>24,:n=>2,:field=>"real",:format=>"array",
                     :author=>"willi",:notes=>"second Line")
end

mktemp() do path, io
    println(io, "%%MatrixMarket Matrix arRAi real GeneraL")
    println(io, "1 2")
    close(io)
    @test_throws DataError MatrixDepot.mmreadheader(path)
end

mktemp() do path, io
    println(io, "%%MatrixMarket Matrix array real GeneraL")
    println(io, "1 2 A")
    close(io)
    @test_throws DataError MatrixDepot.mmreadheader(path)
end

mktemp() do path, io
    println(io, "%%MatrixMarket Matrix coordinate real GeneraL")
    println(io, "1 2")
    close(io)
    @test_throws DataError MatrixDepot.mmreadheader(path)
end

@test MatrixDepot.mmreadheader("") === nothing

# an example loading a txt file
data = mdopen("Pajek/Journals")
@test length(MatrixDepot.metareader(data, "Journals_nodename.txt")) > 100
@test length(MatrixDepot.metareader(data, "Journals_nodename.txt")) > 100 # repeat to use cache

# reading mtx files
io = IOBuffer("""
%%Matrixmarket matrix array real general
2 1
2.1
3.14e0
""")
@test MatrixDepot.mmread(io) == reshape([2.1; 3.14], 2, 1)

io = IOBuffer("""
%%Matrixmarket matrix array real symmetric
2 2
2.1
3.14e0
4
""")
@test MatrixDepot.mmread(io) == [2.1 3.14; 3.14 4.0]

io = IOBuffer("""
%%Matrixmarket matrix array real skew-symmetric
2 2
2.1
3.14e0
""")
@test MatrixDepot.mmread(io) == [0 -2.1; 2.1 0]

# the line containing blanks seemed once to suck when using IOBuffer
io = IOBuffer("%%matrixmarket matrix array complex hermitian\n    \n2 2\n2.1 0\n.314e+1 1\n+3.0 -0.0\n")
@test MatrixDepot.mmread(io) == [2.1 3.14-im; 3.14+im 3]

io = IOBuffer("""
%%Matrixmarket matrix coordinate pattern general
2 2 3
1 2
2 1
2 2
""")
@test MatrixDepot.mmread(io) == [false true; true true]

io = IOBuffer("""
%%Matrixmarket matrix coordinate pattern general
2 2 3
1 2
3 1
2 2
""")
@test_throws DataError MatrixDepot.mmread(io)

io = IOBuffer("""
%%Matricxmarket matrix array real general
2 1
""")
@test_throws Meta.ParseError MatrixDepot.mmread(io)

io = IOBuffer("""
%%Matrixmarket matricx array real general
2 1
""")
@test_throws Meta.ParseError MatrixDepot.mmread(io)

io = IOBuffer("""
%%Matrixmarket matrix arroy real general
2 1
""")
@test_throws Meta.ParseError MatrixDepot.mmread(io)


io = IOBuffer("""
%%Matrixmarket matrix array rehal general
2 1
""")
@test_throws Meta.ParseError MatrixDepot.mmread(io)

io = IOBuffer("""
%%Matrixmarket matrix array real fluffy
2 1
""")
@test_throws Meta.ParseError MatrixDepot.mmread(io)

io = IOBuffer("""
%%Matrixmarket matrix array real general
2 2
1 2
""")
@test_throws BoundsError MatrixDepot.mmread(io)

io = IOBuffer("""
%%Matrixmarket matrix array real symmetric
2 1
1 2
""")
@test_throws DimensionMismatch MatrixDepot.mmread(io)

io = IOBuffer("""
%%Matrixmarket matrix array real symmetric
2 2
1 2
""")
@test_throws BoundsError MatrixDepot.mmread(io)

io = IOBuffer("""
%%Matrixmarket matrix array complex hermitian
2 2
1 2
2 3
""")
@test_throws BoundsError MatrixDepot.mmread(io)

io = IOBuffer("""
%%Matrixmarket matrix array real skew-symmetric
2 2
""")
@test_throws BoundsError MatrixDepot.mmread(io)

