using SynthDiag: add_interferometer!
using SD4SOLPS: preparation
using Test

@testset "interferometer" begin
    eqdsk_file = "g002296.00200"
    sample_dir = "$(@__DIR__)/../sample"
    test_dir = mktempdir()
    ids = preparation(eqdsk_file, sample_dir; filename=test_dir * "/output")
    add_interferometer!("$(@__DIR__)/../src/default_interferometer.json", ids)
    # Just checking if the function runs through for now
    @test true
end
