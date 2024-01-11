using SynthDiag:
    add_interferometer!, add_langmuir_probes!, Noise, add_gas_injection!,
    compute_gas_injection!, instant_gas_model
using OMAS: json2imas, OMAS
using Test
using Printf
using ArgParse: ArgParse
using JSON: JSON
using Plots

function parse_commandline()
    s = ArgParse.ArgParseSettings(; description="Run tests. Default is all tests.")

    ArgParse.add_arg_table!(s,
        ["-i", "--interferometer"],
        Dict(:help => "Test interferometer",
            :action => :store_true),
        ["-g", "--gas_injection"],
        Dict(:help => "Test gas injection",
            :action => :store_true),
        ["-p", "--langmuir_probes"],
        Dict(:help => "Test langmuir probes",
            :action => :store_true),
    )
    args = ArgParse.parse_args(s)
    if !any(values(args)) # If no flags are set, run all tests
        for k ∈ keys(args)
            args[k] = true
        end
    end
    return args
end
args = parse_commandline()

if args["interferometer"]
    @testset "interferometer" begin
        print("json2imas() time: ")
        @time ids =
            json2imas(
                "$(@__DIR__)/../samples/time_dep_edge_profiles_with_equilibrium.json",
            )
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
end

function test_gas_response(config, excitation, plot_title, figname)
    ids = OMAS.dd()
    add_gas_injection!(config, ids)
    ttotal = 5
    nt = Int(ttotal * 1000) + 1
    tstart = 1.0
    tstartind = Int(tstart * 1000) + 1
    tend = 3.0
    tendind = Int(tend * 1000) + 1
    ids.gas_injection.valve[1].voltage.time = collect(LinRange(0, ttotal, nt))
    ids.gas_injection.valve[1].voltage.data = zeros(nt)
    ids.gas_injection.valve[1].voltage.data[tstartind:tendind] .=
        excitation.(ids.gas_injection.valve[1].voltage.time[tstartind:tendind])

    json_data = JSON.parsefile(config)
    valves = json_data["gas_injection"]["valve"]

    ideal_gas_flow_rate = instant_gas_model(
        ids.gas_injection.valve[1].voltage.data,
        valves[1]["dribbling_model"],
    )

    compute_gas_injection!(ids, valves)
    plot(
        ids.gas_injection.valve[1].voltage.time,
        ideal_gas_flow_rate;
        lw=2,
        label="Ideal",
    )
    plot!(
        ids.gas_injection.valve[1].flow_rate.time,
        ids.gas_injection.valve[1].flow_rate.data;
        lw=2,
        label="Dribble model",
    )
    plot!(;
        title=plot_title,
        legend=true,
        xlabel="Time (s)",
        ylabel="Flow rate (Pa m^3/s)",
    )
    return savefig(figname)
end

if args["gas_injection"]
    @testset "gas_injection" begin
        config = "$(@__DIR__)/../src/default_gas_injection.json"
        sine_excitation(t) = 0.3 * cos.(2 * pi * 2 * t) .+ 0.3
        test_gas_response(
            config,
            sine_excitation,
            "Sinusoidal Excitation",
            "$(@__DIR__)/gas_injection_sine.png",
        )
        sine_noise_excitation(t) =
            max(0, 0.2 * cos.(2 * pi * 2 * t) .+ 0.3 + 0.1 * randn())
        test_gas_response(
            config,
            sine_noise_excitation,
            "Noisy Sinusoidal Excitation",
            "$(@__DIR__)/gas_injection_sine_noise.png",
        )
        step_excitation(t) = 0.6
        test_gas_response(
            config,
            step_excitation,
            "Step Excitation",
            "$(@__DIR__)/gas_injection_step.png",
        )
        noise_excitation(t) = max(0.1 * (1 * randn() + 6), 0.0)
        test_gas_response(
            config,
            noise_excitation,
            "Noise Excitation",
            "$(@__DIR__)/gas_injection_noise.png",
        )
        start_stop_excitation(t) = (1.0 < t < 1.5 || 2.0 < t < 2.5) ? 0.6 : 0.0
        test_gas_response(
            config,
            start_stop_excitation,
            "Start-Stop Excitation",
            "$(@__DIR__)/gas_injection_start_stop.png",
        )
        @test true
    end
end

if args["langmuir_probes"]
    @testset "langmuir_probes" begin
        ids =
            json2imas(
                "$(@__DIR__)/../samples/time_dep_edge_profiles_with_equilibrium.json",
            )
        # Assume a 5% noise level in ne values
        ff = 0.0:0.1:1000
        df = ff[2] - ff[1]
        lpf = [f < 10.0 ? 1.0 : 100.0 * f^(-2) for f ∈ ff]
        ne_noise_power = 1e14 * df * lpf
        ne_noise = Noise(ne_noise_power, ff)
        add_langmuir_probes!(
            "$(@__DIR__)/../src/default_langmuir_probes.json",
            ids;
            ne_noise=ne_noise,
        )
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
end
