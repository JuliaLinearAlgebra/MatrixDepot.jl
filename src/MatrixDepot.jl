module MatrixDepot
using Compat # support v3 and v4
using GZip
using MatrixMarket

import Base: show, search

export

# main function
matrixdepot,

@addgroup, @rmgroup

include("ufreader.jl")
include("common.jl")        # main functions
include("higham.jl")        # test matrices
include("regu.jl")          # regularization test problem
include("data.jl")          # matrix data
include("../user/user.jl")  # user defined groups and matrix generators
include("download.jl")      # download data from the UF sparse matrix collection


end # end module
