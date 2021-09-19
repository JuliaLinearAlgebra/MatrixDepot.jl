
const MD_TOGGLE = true

ENV["MATRIXDEPOT_URL_REDIRECT"] = length(VERSION.prerelease) == 0 ? "1" : "1"
const basedir = tempname()
ENV["MATRIXDEPOT_DATA"] = abspath(basedir, "data")
ENV["MATRIXDEPOT_MYDEPOT"] = abspath(basedir, "myMatrixDepot") #Delete when MYDEPOT functionality is deleted

using MatrixDepot
using Test
using LinearAlgebra
using SparseArrays

import MatrixDepot: DataError

function toggle_db()
    MatrixDepot.toggle_remote()
    MatrixDepot.init(ignoredb=true)
    MatrixDepot.update()
end

# that will download the index files if necessary and initialize internal data
MatrixDepot.init()

include("generators.jl")

@testset "MatrixDepot simulate remote matrix tests" begin
    tests = [
            "download",
            "common",
            "number",
            "property",
            ]

    tests = []

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
    toggle_db()
    @test mdopen(sp(1)) != nothing
    mv(xtmpdir, xdir, force=true)
end

@testset "MatrixDepot real remote matrix tests" begin
    toggle_db()
    println("running remote.jl ...")
    include("remote.jl")
    println("finished remote.jl")
end
