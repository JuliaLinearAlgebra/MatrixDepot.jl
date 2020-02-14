
const MD_TOGGLE = true

ENV["MATRIXDEPOT_URL_REDIRECT"] = length(VERSION.prerelease) == 0 ? "1" : "1"
const basedir = tempname()
ENV["MATRIXDEPOT_DATA"] = abspath(basedir, "data")
ENV["MATRIXDEPOT_MYDEPOT"] = abspath(basedir, "myMatrixDepot")

using MatrixDepot
using Test
using LinearAlgebra
using SparseArrays

import MatrixDepot: DataError

# include("clean.jl")
# data_save = save_target(data_dir)
# user_save = save_target(user_dir)

# that will download the index files if necessary and initialize internal data
MatrixDepot.init()

include("generators.jl")

@testset "MatrixDepot remote matrix tests" begin
    tests = [
            "download",
            "common",
            "number",
            "property",
            ]

    @testset "$t" for t in tests
        tp = joinpath(@__DIR__(), "$(t).jl")
        println("running $(tp) ...")
        include(tp)
        println("finished $(tp)")
    end

    xdir = MatrixDepot.data_dir()
    xtmpdir = string(xdir, ".tmp")
    println("mv $xdir $xtmpdir")
    mv(xdir, xtmpdir, force=true)
    MatrixDepot.toggle_remote()
    MatrixDepot.init(ignoredb=true)
    MatrixDepot.update()
    @test mdopen(sp(1)) != nothing
    mv(xtmpdir, xdir, force=true)
end

