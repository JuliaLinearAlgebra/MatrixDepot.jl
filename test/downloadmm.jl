
@testset "mdalternative interface" begin
    mpat = "Harwell-Boeing/psadmit/1138_bus"
    spat = "HB/1138_bus"
    @test mdlist(mm(spat)) == [mpat]
    @test mdlist(sp(mpat)) == [spat]
    @test mdlist(mm(mpat)) == [mpat]
    @test mdlist(sp(spat)) == [spat]
    @test mdlist(sp("vand")) == String[]
    @test mdlist(mm("Sandia/adder_dcop_17")) == String[]

    s1 = mdlist(sp(mm(:)))
    @test length(s1) == 455
    m1 = mdlist(mm(s1))
    @test length(m1) == 455
    @test m1 ⊆ mdlist(mm(:))
    
    m2 = mdlist(mm(sp(:)))
    @test length(m2) == 455
    s2 = mdlist(sp(m2))
    @test length(s2) == 455
    @test s2 ⊆ mdlist(sp(:))

    @test m1 == m2
end
