
ENV["MATRIXDEPOT_URL_REDIRECT"] = length(VERSION.prerelease) == 0 ? "1" : "1"
const basedir = tempname() # tests start with a clean data directory
ENV["MATRIXDEPOT_DATA"] = abspath(basedir, "data")

using MatrixDepot
using Test
using LinearAlgebra
using SparseArrays

import MatrixDepot: DataError

# that will download the index files if necessary and initialize internal data
MatrixDepot.init()

@testset "MatrixDepot.jl" begin

include("generators.jl")
include("include_generator.jl")

@testset "MatrixDepot simulate remote matrix tests" begin
    tests = [
            "download",
            "common",
            "number",
            "property",
            "downloadmm",
            ]

    @testset "$t" for t in tests
        tp = joinpath(@__DIR__(), "$(t).jl")
        println("running $(tp) ...")
        include(tp)
        println("finished $(tp)")
    end
end

@testset "MatrixDepot real remote matrix tests" begin
    println("running remote.jl ...")
    include("remote.jl")
    println("finished remote.jl")
end

end # testset MatrixDepot.jl
