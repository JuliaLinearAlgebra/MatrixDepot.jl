dirdata = joinpath(dirname(@__FILE__),"..","data")
if isdir(string(dirdata, '/', "uf"))
    rm(string(dirdata, '/', "uf"), recursive = true)
end
if isdir(string(dirdata, '/', "mm"))
    rm(string(dirdata, '/', "mm"), recursive = true)
end

#download
# MatrixDepot.update() # already done when initializing directory
# uf
matrixdepot("1138_bus", :get)
matrixdepot("HB/1138_bus", :get)
matrixdepot("Pajek/Journals", :get)
# Matrix Markt
matrixdepot("Harwell-Boeing/smtape/bp___200", :get)
 
# read data
A = matrixdepot("HB/1138_bus", :r)
matrixdepot("HB/1138_bus")
matrixdepot("1138_bus", :s) #search
#B = matrixdepot("Harwell-Boeing/psadmit/662_bus", :read)
#matrixdepot("Harwell-Boeing/psadmit/662_bus")
r = matrixdepot("Pajek/Journals", :read, meta = true)
display(r["Journals"])
# MatrixDepot.update() # already done

matrixdepot("Bates/*", :get)
B = matrixdepot("Bates/Chem97Zt", :r)
matrixdepot("Bates/Chem97Zt")
matrixdepot()
@test_throws ArgumentError matrixdepot("HB/662_bus", :k)

matrixdepot("epb0", :get)

# matrix market
matrixdepot("Harwell-Boeing/lanpro/nos5", :get)
matrixdepot()

@testset "load and read metadata" begin
# an example with rhs and solution
data = mdopen("DRIVCAV/cavity14")
@test size(matrix(data)) == (2597, 2597)
@test size(rhs(data), 1) == 2597
@test size(solution(data), 1) == 2597

# read a format array file
@test MatrixDepot.fileinfo(abspath(MatrixDepot.matrixfile(data), "..", string("cavity14_b.mtx"))) != nothing

# read from a pipeline
@test open(`echo -e '%%matrixmarket matrix array real general\n1 1\n2.5'`) do io
    MatrixDepot.mmread(io) == reshape([2.5], 1, 1)
end

# an example loading a txt file
data = mdopen("Pajek/Journals")
@test length(metareader(data, "Journals_nodename.txt")) > 100

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

io = IOBuffer("""
%%Matrixmarket matrix array complex herMitian
2 2
2.1 0
3.14e0 1.0
3.0 0
""")
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
@test_throws ArgumentError MatrixDepot.mmread(io)

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




end
# rm data

#=
rm(string(dirdata, '/', "uf_matrices.html"))
rm(string(dirdata, '/', "mm_matrices.html"))
rm(string(dirdata, '/', "uf"), recursive = true)
rm(string(dirdata, '/', "mm"), recursive = true)
=#

