module MatrixDepot
using Compat # support v3 and v4
using GZip
using MatrixMarket

export

# main function
matrixdepot,

@addproperty, @rmproperty

include("ufreader.jl")
include("common.jl")        # main functions
include("higham.jl")        # test matrices
include("data.jl")
include("user.jl")          # user defined properties
include("download.jl")      # download data from the UF sparse matrix collection


end # end module
