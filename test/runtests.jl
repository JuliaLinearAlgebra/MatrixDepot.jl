using MatrixDepot
using Base.Test

include("test_cauchy.jl")
include("test_circul.jl")
include("test_hadamard.jl")
include("test_hilb.jl")
include("test_dingdong.jl")
include("test_frank.jl")
include("test_invhilb.jl")

println("Success in all tests.")
