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
@test loadsvd("HB/1138_bus") == 0 # already loaded
@test (@test_logs (:info, r"downloading") loadsvd("HB/494_bus")) == 0 # not available
@test loadinfo("HB/494_bus") == 0
@test MatrixDepot.addsvd!(MatrixDepot.mdata("HB/494_bus")) == 0
@test loadsvd("hilb") == loadinfo("hilb") == 0
@test (@test_logs (:warn, r"could not load hilb") MatrixDepot._load(MatrixDepot.MATRIX_DB, error, "hilb")) == 0

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

data = MatrixDepot.mdata("HB/1138_bus")
@test (@test_logs (:warn, r"quadratic") MatrixDepot.extremnnz(data, 2, 3, 4)) == (5, 8)
kind = data.kind
data.header.kind = "other kind"
@test_logs (:warn, r"header kind") match_mode=:any MatrixDepot.addmetadata!(data)
data.header.kind = kind
@test MatrixDepot.addmetadata!(MatrixDepot.mdata("hilb")) == 0
@test !MatrixDepot.issvdok(MatrixDepot.mdata("hilb"))
@test MatrixDepot._vector(2.0) == [2.0]
@test MatrixDepot.svdstatistics([1.0, 1e-20, 0.0]).rank == 1
@test occursin(" ... ", sprint(show, mdopen("*/bfly").data))
@test sprint(show, MatrixDepot.MATRIX_DB) == "MatrixDatabase($(length(MatrixDepot.MATRIX_DB.data)))"

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
mktemp() do path, io
    close(io)
    @test MatrixDepot.mmreadheader(path) === nothing
    @test_throws DataError MatrixDepot.mmread(path)
end
mktemp() do path, io # header of a second matrix is used when the first has no nz
    println(io, "%%MatrixMarket matrix array real general\n1 1\n1.0")
    println(io, "%%MatrixMarket matrix coordinate real general\n2 2 1\n1 1 1.0")
    close(io)
    @test MatrixDepot.mmreadheader(path)[:nz] == 1
end

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

A = sparse([1, 2, 2], [1, 1, 2], [1.0, 2im, 3.0])
@test MatrixDepot.mtranspose(A) == transpose(A)
@test MatrixDepot.madjoint(A) == A'
@test MatrixDepot.mtranspose(1:2) == transpose(1:2)
@test MatrixDepot.mksparse(2, 2, [1, 2], [2, 1], [1.0, 2.0]) == [0 1; 2 0]

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

# last entry is not followed by a newline
mktempdir() do dir
    for (text, A) in (
        ("coordinate real general\n2 2 1\n1 1 1.5", sparse([1], [1], [1.5], 2, 2)),
        ("coordinate integer general\n2 2 1\n1 1 7", sparse([1], [1], [7], 2, 2)),
        ("coordinate pattern general\n2 2 1\n2 1", sparse([2], [1], [true], 2, 2)),
        ("coordinate complex general\n2 2 1\n1 2 1 -2", sparse([1], [2], [1-2im], 2, 2)),
        ("coordinate real general\n2 2 2\n1 1 1.5\n2 2 -2.5e1", sparse([1, 2], [1, 2], [1.5, -25], 2, 2)),
        ("array real general\n2 1\n1.5\n2.5", reshape([1.5, 2.5], 2, 1)),
        ("array integer general\n1 2\n3\n-4", [3 -4]),
        ("array complex general\n1 1\n1.5 2", fill(1.5+2im, 1, 1)),
        )
        path = joinpath(dir, "nonewline.mtx")
        write(path, "%%MatrixMarket matrix " * text)
        @test MatrixDepot.mmread(path) == A
        @test MatrixDepot.mmread(IOBuffer(read(path))) == A
    end
end

# file ends before the size line - used to loop forever, so run in a subprocess with time limit
mktempdir() do dir
    files = [joinpath(dir, "header.mtx"), joinpath(dir, "comments.mtx")]
    write(files[1], "%%MatrixMarket matrix coordinate real general\n")
    write(files[2], "%%MatrixMarket matrix array real general\n% comment\n\n%\n")
    code = """
    using MatrixDepot
    for f in ARGS, g in (MatrixDepot.mmread, MatrixDepot.mmreadheader)
        try g(f); println("no exception") catch ex; println(nameof(typeof(ex))) end
    end
    """
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $code $files`
    out = Pipe()
    proc = run(pipeline(cmd, stdout=out), wait=false)
    close(out.in)
    output = @async read(out, String)
    timedwait(() -> process_exited(proc), 300.0) == :ok || kill(proc)
    @test success(proc)
    @test split(fetch(output)) == repeat(["ParseError", "DataError"], length(files))
end

# corrupted database file is recreated
db = MatrixDepot.MATRIX_DB
write(MatrixDepot.dbpath(db), "garbage")
@test_logs (:warn, r"recreating") match_mode=:any MatrixDepot.downloadindices(db)
