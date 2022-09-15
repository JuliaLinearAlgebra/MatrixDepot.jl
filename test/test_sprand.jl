
@testset "corner cases sprandsvd" begin
    A = matrixdepot("sprandsvd", 0, 0.1, 10.0)
    @test A isa SparseMatrixCSC
    @test A == spzeros(0, 0, 0.0)
end
