
import MatrixDepot: namemm2ss, namess2mm, Pattern


function namex(x::AbstractString)
    n = length(split(x, '/'))
    y = n == 2 ? namess2mm(x) : n == 3 ? namemm2ss(x) : String[]
    x != y ? y : String[]
end

"""
    mdalternative(p::Pattern)

Like `mdlist`, but return list of problem names from alternative source.

Most of the problems of the `Matrixmarket` collection are also available in the `Suite Sparse` collection,
sometimes with modified name, sometimes also the content is slightly different.

If no alternative is present, an empty list is returned.
"""
function mdalternative(p::Pattern)
    vcat(String[], (mdlist(namex(x)) for x in mdlist(p))...)
end


@testset "provisional mdalternative interface" begin
    s1 = mdalternative(mm(:))
    @test length(s1) == 455
    m1 = mdalternative(s1)
    @test m1 ⊆ mdlist(mm(:))
    
    m2 = mdalternative(sp(:))
    @test length(m2) == 455
    s2 = mdalternative(m2)
    @test s2 ⊆ mdlist(sp(:))
end
