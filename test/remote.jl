#download from remote site
import MatrixDepot: URL_REDIRECT, load, loadinfo, loadsvd

URL_REDIRECT[] = 0
# sp
pattern = "*/494_bus"
# download.jl marked the SVD data of HB/494_bus as not available; start clean
rm(MatrixDepot.localdir(MatrixDepot.mdata("HB/494_bus")), force=true, recursive=true)
@test loadinfo(pattern) == 1 # load only header
println("downloading svd for $pattern")
@test loadsvd(pattern) == 1 # load svd extension data
@test load(pattern) == 1 # loaded full data

# Matrix Market
pattern = "*/*/494_bus" # which has not been touched in other tests
@test loadinfo(pattern) == 1 # load only header
@test load(pattern) == 1 # loaded full data
