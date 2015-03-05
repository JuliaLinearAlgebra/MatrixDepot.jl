module MatrixDepot
using Compat # support v3 and v4 

export matrixdepot, @addproperty, @rmproperty

include("matrixdepot.jl")    # main functions 
include("higham.jl")      # test matrices
include("user.jl")        # user defined properties


end # end module
