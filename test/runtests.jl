
const MD_TOGGLE = true

using MatrixDepot
using Test
using LinearAlgebra
using SparseArrays

import MatrixDepot: DataError

macro inc(a)
    :(@testset $a begin let n, p, i; include($a) end end)
end

include("clean.jl")
data_save = save_target(data_dir)
user_save = save_target(user_dir)

# that will download the index files if necessary and initialize internal data
MatrixDepot.init()
@testset "MatrixDepot generator tests" begin

    @inc("test_magic.jl")
    @inc("test_cauchy.jl")
    @inc("test_circul.jl")
    @inc("test_hadamard.jl")
    @inc("test_hilb.jl")
    @inc("test_dingdong.jl")
    @inc("test_frank.jl")
    @inc("test_invhilb.jl")
    @inc("test_forsythe.jl")
    @inc("test_grcar.jl")
    @inc("test_triw.jl")
    @inc("test_moler.jl")
    @inc("test_pascal.jl")
    @inc("test_kahan.jl")
    @inc("test_pei.jl")
    @inc("test_vand.jl")
    @inc("test_invol.jl")
    @inc("test_chebspec.jl")
    @inc("test_lotkin.jl")
    @inc("test_clement.jl")
    @inc("test_fiedler.jl")
    @inc("test_minij.jl")
    @inc("test_binomial.jl")
    @inc("test_tridiag.jl")
    @inc("test_lehmer.jl")
    @inc("test_parter.jl")
    @inc("test_chow.jl")
    @inc("test_randcorr.jl")
    @inc("test_poisson.jl")
    @inc("test_neumann.jl")
    @inc("test_rosser.jl")
    @inc("test_sampling.jl")
    @inc("test_wilkinson.jl")
    @inc("test_rando.jl")
    @inc("test_randsvd.jl")
    @inc("test_rohess.jl")
    @inc("test_kms.jl")
    @inc("test_wathen.jl")
    @inc("test_oscillate.jl")
    @inc("test_toeplitz.jl")
    @inc("test_hankel.jl")
    @inc("test_golub.jl")
    @inc("test_companion.jl")
    @inc("test_erdrey.jl")
    @inc("test_gilbert.jl")
    @inc("test_smallworld.jl")
end

@testset "MatrixDepot remote matrix tests" begin
    tests = [
            "regu",
            "include_generator",
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
        @test mdopen(uf(1)) != nothing
        mv(xtmpdir, xdir, force=true)
    finally
        revert_target(user_save, user_dir)
        revert_target(data_save, data_dir)
    end
end

println("Success in all tests.")
