module MatrixDepot
using Compat # support v3 and v4 
using MAT

export 

# main function
matrixdepot, 

@addproperty, @rmproperty,

# download
downloadsparse, updatesparse

include("mmreader.jl")      # Matrix Market format reader  
include("matrixdepot.jl")   # main functions 
include("higham.jl")        # test matrices
include("user.jl")          # user defined properties
include("download.jl")      # download data from the UF sparse matrix collection

end # end module
