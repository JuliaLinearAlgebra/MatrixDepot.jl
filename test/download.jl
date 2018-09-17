#download
# verify that no data are downloaded (empty directory)
@test list(isloaded) == []

# uf
@test load("**/1138_bus") == 2 # loaded both versions (uf and mm)
@test load("HB/1138_bus") == 0 # count actually loaded
@test load("Pajek/Journals") == 1
@test load("wing") == 0 # do not try to load local matrix
# Matrix Markt
@test load("Harwell-Boeing/smtape/bp___200") == 1

@test list(isloaded) == ["HB/1138_bus", "Harwell-Boeing/psadmit/1138_bus",
                        "Harwell-Boeing/smtape/bp___200", "Pajek/Journals"]

# read data
@test matrixdepot("HB/1138_bus") != nothing
@test list("**/1138_bus") == ["HB/1138_bus", "Harwell-Boeing/psadmit/1138_bus"]
@test string(mdinfo("Harwell-Boeing/psadmit/662_bus")) ==
        "# Harwell-Boeing/psadmit/662_bus\n\n" *
        "###### MatrixMarket matrix coordinate real symmetric\n\n662 662 1568\n"

@test MatrixDepot.loadinfo(MatrixDepot.mdata("HB/1138_bus")) == 0
@test MatrixDepot.loadinfo(MatrixDepot.mdata("Bates/Chem97Zt")) == 1

# metatdata access with old interface
@test metareader(mdopen("Pajek/Journals"), "Journals.mtx") !== nothing
@test_throws DataError metareader(mdopen("*/1138_bus"), "1138_bus_b")
@test_throws DataError metareader(mdopen("that is nothing"))

@test load("Bates/C*") == 2
@test mdinfo("Bates/Chem97Zt") != nothing
@test matrixdepot("Bates/Chem97Zt") != nothing
@test "Bates/Chem97Zt" in list(isloaded)
@test length(list(isloaded)) == 6

@test load("**/epb0") == 1

# matrix market
@test load("Harwell-Boeing/lanpro/nos5") == 1
@test length(list(isloaded)) == 8
@test mdopen("Bai/dwg961b") !== nothing
mdesc = mdopen("Bai/dwg961b")
@test iscomplex(mdesc.data)
@test issymmetric(mdesc.data)
@test !ishermitian(mdesc.data)
@test !ispattern(mdesc.data)
@test metadata(mdesc) == ["dwg961b.mtx"]
@test metareader(mdesc, "dwg961b.mtx") == matrixdepot("Bai/dwg961b")

# an example with rhs and solution
mdesc = mdopen("DRIVCAV/cavity14")
A1 = mdesc.A
@test size(A1) == (2597, 2597)
A2 = mdesc.A
@test A1 === A2 # cache should be used
@test size(mdesc.b, 1) == 2597
@test size(mdesc.x, 1) == 2597
@test_throws DataError metareader(mdesc, "invlid")
@test_throws DataError mdesc.y

# read a format array file
@test MatrixDepot.mmreadheader(abspath(MatrixDepot.matrixfile(mdesc.data),
                                       "..", string("cavity14_b.mtx"))) != nothing

mdesc = mdopen("blur", 10, false)
@test mdesc.A isa AbstractMatrix
@test mdesc.b isa AbstractVector
@test mdesc.x isa AbstractVector
@test_throws DataError metareader(mdesc, "invlid")
@test_throws DataError mdesc.y

@test_throws DataError mdopen(uf(1)).b # no rhs data for this example

mdesc = mdopen("*/bfly")
@test metareader(mdesc, "bfly_Gname_01.txt") == "BFLY3\n"
@test metareader(mdesc, "Gname_01.txt") == "BFLY3\n"
@test size(mdesc.G_06) == (2048, 2048)
@test mdesc.G_06 === metareader(mdesc, "G_06")
@test mdesc.G_06 === metareader(mdesc, "bfly_G_06.mtx")
fn = joinpath(dirname(MatrixDepot.matrixfile(mdesc.data)), "bfly_Gname_01.txt")
@test_throws DataError MatrixDepot.mmreadheader(fn)

# read from a pipeline
@test open(`printf '%%%%matrixmarket matrix array real general\n1 1\n2.5\n'`) do io
    MatrixDepot.mmread(io) == reshape([2.5], 1, 1)
end

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
@test length(metareader(data, "Journals_nodename.txt")) > 100
@test length(metareader(data, "Journals_nodename.txt")) > 100 # repeat to use cache

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

# the line containing blanks seemed to suck when using IOBuffer
@test open(`printf '%%%%matrixmarket matrix array complex hermitian\n    \n2 2\n2.1 0\n.314e+1 1\n+3.0 -0.0\n'`) do io
    MatrixDepot.mmread(io) == [2.1 3.14-im; 3.14+im 3]
end

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

