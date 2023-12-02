using SynthDiag: add_interferometer!, add_langmuir_probes!
using OMAS: json2imas
using Test
using Printf

@testset "interferometer" begin
    ids =
        json2imas("$(@__DIR__)/../samples/time_dep_edge_profiles_with_equilibrium.json")
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

@testset "langmuir_probes" begin
    ids =
        json2imas("$(@__DIR__)/../samples/time_dep_edge_profiles_with_equilibrium.json")
    add_langmuir_probes!("$(@__DIR__)/../src/default_langmuir_probes.json", ids)
    # Just checking if the function runs through for now
    for lp ∈ ids.langmuir_probes.embedded
        println()
        println("-"^49)
        println("-"^49)
        println("Probe: $(lp.name)")
        println("-"^49)
        @printf("|%15s|%15s|%15s|\n", "time", "n_e", "t_e")
        println("-"^49)
        for ii ∈ eachindex(lp.time)
            @printf(
                "|%15.3e|%15.3e|%15.3e|\n",
                lp.time[ii],
                lp.n_e.data[ii],
                lp.t_e.data[ii]
            )
        end
        println("-"^49)
    end
    @test true
end