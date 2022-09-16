
@testset "SparseMatrixD" begin
    A = sprand(10, 20, 0.2)
    D = MatrixDepot.SparseMatrixD(A)
    io = IOBuffer()
    show(io, MIME("text/plain"), D)
    @test size(D) == (10, 20)
    @test A == D
end

@testset "corner cases sprandsvd ($m, $n)" for (m, n) in [(0, 0); (0, 10); (20, 0)]
    A = matrixdepot("sprandsvd", m, n, 0.1)
    @test A isa SparseMatrixCSC
    @test A == spzeros(m, n)
end
@testset "corner cases sprandsvd ($m, $n, 10.0, $m )" for (m, n) in [(1, 1);], mode in 1:5
    A = matrixdepot("sprandsvd", m, n, 0.1, 10.0)
    @test A isa SparseMatrixCSC
    @test size(A) == (m, n)
    @test A[1,1] == 10.0
end
@testset "corner cases sprandsym ($n)" for n in [0;]
    A = matrixdepot("sprandsvd", n, 0.1)
    @test size(A) == (n, n)
    A = matrixdepot("sprandsym", n, 0.1, 10.0)
    @test A isa SparseMatrixCSC
    @test A == spzeros(n, n)
end
@testset "corner cases sprandsym ($n, 10.0, $m )" for n in [1;], mode in 1:5
    A = matrixdepot("sprandsym", n, 0.1, 10.0, mode)
    @test A isa SparseMatrixCSC
    @test size(A) == (n, n)
    @test A[1,1] == 10.0
end

@testset "$fname argument errors $fname" for fname in ["sprandsvd", "sprandsym"]
    n = 10
    @test_throws ArgumentError matrixdepot(fname, n, 1.1, 10.0)
    @test_throws ArgumentError matrixdepot(fname, n, 0.1, 0.9)
    @test_throws ArgumentError matrixdepot(fname, n, 0.1, 10.0, 0)
    @test_throws MethodError matrixdepot(fname, 1, 1, 1, 1, 1, 1, 1)
end

@testset "sprandsvd ($m, $n, $p)" for (m, n) in [(10, 10); (10, 12); (20, 10)], p in [0.0; 0.2; 1.0]
    A = matrixdepot("sprandsvd", m, n, p, 10.0, 4)
    @test svdvals(Matrix(A)) ≈ 1.0:-0.1:0.1
end

@testset "sprandsvd ($m, $n, $p, $s)" for (m, n) in [(10, 5)], p in [0.5], s in [[1,2], [1,2,3,4,5]]
    A = matrixdepot("sprandsvd", m, n, p, s)
    @test A isa SparseMatrixCSC
    @test size(A) == (m, n)
    @test nnz(A) >= Int(floor(p * m * n))
    es = sort!([s; zeros(min(m, n)-length(s))], rev=true)
    @test svdvals(Matrix(A)) ≈ es
end
@testset "sprandsvd (ComplexF64, $m, $n, $p, $s)" for (m, n) in [(10, 5)], p in [0.5], s in [[1,2], [1,2,3,4,5]]
    A = matrixdepot("sprandsvd", m, n, p, s)
    @test A isa SparseMatrixCSC
    @test size(A) == (m, n)
    @test nnz(A) >= Int(floor(p * m * n))
    es = sort!([s; zeros(min(m, n)-length(s))], rev=true)
    @test svdvals(Matrix(A)) ≈ es
end

@testset "sprandsym ($n, $p, $s)" for n in [5], p in [0.5], s in [[-2,2], [-1,-2,0,4,5]]
    A = matrixdepot("sprandsym", n, p, s)
    @test A isa SparseMatrixCSC
    @test size(A) == (n, n)
    @test nnz(A) >= Int(floor(p * m * n))
    @test issymmetric(A)
    es = sort!([s; zeros(min(m, n)-length(s))], rev=false)
    @test eigvals(Matrix(A)) ≈ es
end
@testset "sprandsym ($n, $p, $s)" for n in [5], p in [0.5], s in [[-2,2], [-1,-2,0,4,5]]
    A = matrixdepot("sprandsym", ComplexF64, n, p, s)
    @test A isa SparseMatrixCSC
    @test size(A) == (n, n)
    @test nnz(A) >= Int(floor(p * m * n))
    @test ishermitian(A)
    es = sort!([s; zeros(min(m, n)-length(s))], rev=false)
    @test eigvals(Matrix(A)) ≈ es
end

println("'sprand...' passed test...")
