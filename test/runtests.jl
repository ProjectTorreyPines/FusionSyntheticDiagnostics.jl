using SynthDiag: add_interferometer!
using OMAS: json2imas
using Test
using Printf

@testset "interferometer" begin
    # eqdsk_file = "g002296.00200"
    # sample_dir = "$(@__DIR__)/../sample"
    # test_dir = mktempdir()
    # ids = preparation(eqdsk_file, sample_dir; filename=test_dir * "/output")
    ids =
        json2imas("$(@__DIR__)/../sample/time_dep_edge_profiles_with_equilibrium.json")
    add_interferometer!("$(@__DIR__)/../src/default_interferometer.json", ids)
    # Just checking if the function runs through for now
    for ch ∈ ids.interferometer.channel
        println()
        println("-"^49)
        println("-"^49)
        println("Channel $(ch.name)")
        println("-"^49)
        @printf("|%15s|%15s|%15s|\n", "time", "int_n_e", "n_e_average")
        println("-"^49)
        for ii ∈ eachindex(ch.n_e_line.data)
            @printf(
                "|%15.3e|%15.3e|%15.3e|\n",
                ch.n_e_line.time[ii],
                ch.n_e_line.data[ii],
                ch.n_e_line_average.data[ii]
            )
        end
        println("-"^49)
    end
    @test true
end
