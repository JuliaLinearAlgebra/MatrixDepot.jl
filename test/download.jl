#download
# verify that no data are downloaded (empty directory)
@test list(:loaded) == []

# uf
@test matrixdepot("1138_bus", :get) == 2 # loaded both versions (uf and mm)
@test matrixdepot("HB/1138_bus", :get) == 1
@test matrixdepot("Pajek/Journals", :get) == 1
# Matrix Markt
@test matrixdepot("Harwell-Boeing/smtape/bp___200", :get) == 1

@test list(:loaded) == ["HB/1138_bus", "Harwell-Boeing/psadmit/1138_bus",
                        "Harwell-Boeing/smtape/bp___200", "Pajek/Journals"]

# read data
@test matrixdepot("HB/1138_bus", :r) != nothing
@test matrixdepot("HB/1138_bus") != nothing
@test matrixdepot("1138_bus", :s) == ["HB/1138_bus", "Harwell-Boeing/psadmit/1138_bus"]
@test_throws DataError matrixdepot("Harwell-Boeing/psadmit/662_bus", :read)
@test string(matrixdepot("Harwell-Boeing/psadmit/662_bus")) ==
                         "# Harwell-Boeing/psadmit/662_bus\n\nno info available\n"

# metatdata access with old interface
@test matrixdepot("Pajek/Journals", :read, meta = true)["Journals"] !== nothing
@test_throws DataError matrixdepot("1138_bus", :read, meta = true)
@test_throws DataError matrixdepot("that is nothing", :read, meta = true)

@test matrixdepot("Bates/C*", :get) == 2
@test matrixdepot("Bates/Chem97Zt") != nothing
@test matrixdepot("Bates/Chem97Zt", :r) != nothing
@test "Bates/Chem97Zt" in list(:loaded)
@test length(list(:loaded)) == 6

@test_throws ArgumentError matrixdepot("HB/662_bus", :k)

@test matrixdepot("epb0", :get) == 1

# matrix market
@test matrixdepot("Harwell-Boeing/lanpro/nos5", :get) == 1
@test length(list(:loaded)) == 8
@test mdopen("Bai/dwg961b") !== nothing
data = mdopen("Bai/dwg961b")
@test iscomplex(data)
@test issymmetric(data)
@test !ishermitian(data)
@test !ispattern(data)
@test metadata(data) == ["dwg961b.mtx"]

# an example with rhs and solution
data = mdopen("DRIVCAV/cavity14"; cache=true)
@test size(matrix(data)) == (2597, 2597)
@test size(matrix(data)) == (2597, 2597) # cache should be used (coverage)
@test size(rhs(data), 1) == 2597
@test size(solution(data), 1) == 2597
@test mdclose(data) === data
@test_throws DataError metareader(data, "invlid")

@test_throws DataError rhs(uf(1)) # no rhs data for this example

# read a format array file
@test MatrixDepot.fileinfo(abspath(MatrixDepot.matrixfile(data), "..", string("cavity14_b.mtx"))) != nothing

data = mdopen("*/bfly")
@test metareader(data, "bfly_Gname_01.txt") == "BFLY3\n"
fn = joinpath(dirname(MatrixDepot.matrixfile(data)), "bfly_Gname_01.txt")
@test_throws DataError MatrixDepot.fileinfo(fn)

# read from a pipeline
@test open(`printf '%%%%matrixmarket matrix array real general\n1 1\n2.5\n'`) do io
    MatrixDepot.mmread(io) == reshape([2.5], 1, 1)
end

# read from temporary file
mktemp() do path, io
    println(io, "%%MatrixMarket Matrix arRAy real GeneraL")
    println(io, "% first line")
    println(io)
    println(io, "%second line")
    println(io); println(io, "    ")
    println(io, "24 2")
    close(io)
    @test MatrixDepot.ufinfo(path) ==
    "%%MatrixMarket Matrix arRAy real GeneraL\n% first line\n\n%second line\n\n    \n24 2\n"
    @test MatrixDepot.fileinfo(path) == (24, 2, 0, "matrix", "array", "real", "general")
end

mktemp() do path, io
    println(io, "%%MatrixMarket Matrix arRAi real GeneraL")
    println(io, "1 2")
    close(io)
    @test_throws DataError MatrixDepot.fileinfo(path)
end

mktemp() do path, io
    println(io, "%%MatrixMarket Matrix array real GeneraL")
    println(io, "1 2 A")
    close(io)
    @test_throws DataError MatrixDepot.fileinfo(path)
end

mktemp() do path, io
    println(io, "%%MatrixMarket Matrix coordinate real GeneraL")
    println(io, "1 2")
    close(io)
    @test_throws DataError MatrixDepot.fileinfo(path)
end

@test MatrixDepot.fileinfo("") === nothing

# an example loading a txt file
data = mdopen("Pajek/Journals"; cache=true)
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

