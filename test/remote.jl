
#download from remote site
import MatrixDepot: URL_REDIRECT, load, loadinfo, loadsvd

# test only for julia nightly
if length(VERSION.prerelease) == 0
    println("remote tests not executed")
else
    URL_REDIRECT[] = 0
    # sp
    pattern = "*/494_bus" # which has not been touched in other tests
    @test loadinfo(pattern) == 1 # load only header
    println("downloading svd for $pattern")
    @test loadsvd(pattern) == 1 # load svd extension data
    @test load(pattern) == 1 # loaded full data

    # Matrix Market
    pattern = "*/*/494_bus" # which has not been touched in other tests
    @test loadinfo(pattern) == 1 # load only header
    @test load(pattern) == 1 # loaded full data
end
