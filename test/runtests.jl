
const MD_TOGGLE = true

using MatrixDepot
using Test
using LinearAlgebra
using SparseArrays

import MatrixDepot: DataError

MatrixDepot.URL_REDIRECT[] = true

include("clean.jl")
data_save = save_target(data_dir)
user_save = save_target(user_dir)

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

    try
        @testset "$t" for t in tests
            tp = joinpath(@__DIR__(), "$(t).jl")
            println("running $(tp) ...")
            include(tp)
            println("finished $(tp)")
        end

        xdir = MatrixDepot.DATA_DIR
        xtmpdir = string(xdir, ".tmp")
        println("mv $xdir $xtmpdir")
        mv(xdir, xtmpdir, force=true)
        MatrixDepot.toggle_remote()
        MatrixDepot.init(ignoredb=true)
        MatrixDepot.update()
        @test mdopen(sp(1)) != nothing
        mv(xtmpdir, xdir, force=true)
    finally
        revert_target(user_save, user_dir)
        revert_target(data_save, data_dir)
    end
end

println("Success in all tests.")
