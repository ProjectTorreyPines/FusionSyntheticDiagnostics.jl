using SynthDiag: IMAS, add_interferometer!, add_langmuir_probes!, add_gas_injection!,
    compute_gas_injection!, get_gas_injection_response, Noise, OverwriteAttemptError,
    langmuir_probe_current, magic_nesep,
    calc_loss_power, calc_conducted_loss_power, calc_q_cyl, calc_heat_flux_width,
    find_OMP_RZ, read_B_theta_OMP, summarize_flux_surfaces!, read_B_theta_OMP_no_ggd
using IMAS: json2imas, gradient
using Test
using Printf
using Plots
using ArgParse: ArgParse
using DelimitedFiles: readdlm
using Interpolations: linear_interpolation

function parse_commandline()
    # Define newARGS = ["--yourflag"] to run only tests on your flags when including runtests.jl
    localARGS = @isdefined(newARGS) ? newARGS : ARGS  # Thanks https://stackoverflow.com/a/44978474/6605826
    s = ArgParse.ArgParseSettings(; description="Run tests. Default is all tests.")

    ArgParse.add_arg_table!(s,
        ["--interferometer"],
        Dict(:help => "Test only interferometer",
            :action => :store_true),
        ["--langmuir_probes"],
        Dict(:help => "Test only langmuir probes",
            :action => :store_true),
        ["--gas_injection"],
        Dict(:help => "Test only gas injection",
            :action => :store_true),
        ["--magic"],
        Dict(:help => "Test only magic diagnostics",
            :action => :store_true),
        ["--derived"],
        Dict(:help => "Test only derived quantities and calculations",
            :action => :store_true),
    )
    args = ArgParse.parse_args(localARGS, s)
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
        ids =
            json2imas(
                "$(@__DIR__)/../samples/time_dep_edge_profiles_with_equilibrium.json",
            )
        add_interferometer!(
            "$(@__DIR__)/../src/default_interferometer.json",
            ids;
            n_e_gsi=-5,
        )
        # Just checking if the function runs through for now
        for ch ∈ ids.interferometer.channel
            println()
            println("-"^49)
            println("-"^49)
            println("Channel $(ch.name)")
            println("-"^49)
            @printf("|%15s|%15s|%15s|\n", "time", "int_n_e", "n_e_average")
            println("-"^49)
            for ii ∈ 1:20:length(ch.n_e_line.time)
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
        # Try overwriting the interferometer data
        try
            add_interferometer!(
                "$(@__DIR__)/../samples/test_interferometer_same_names.json",
                ids; n_e_gsi=-5,
            )
        catch err
            @test isa(err, OverwriteAttemptError)
        end
        # Try overwriting the interferometer data with the overwrite flag
        add_interferometer!(
            "$(@__DIR__)/../samples/test_interferometer_same_names.json",
            ids;
            overwrite=true, n_e_gsi=-5,
        )
        @test ids.interferometer.channel[1].name == "V1"
        @test ids.interferometer.channel[1].line_of_sight.first_point.r == 5.2
        @test ids.interferometer.channel[2].name == "V2"
        @test ids.interferometer.channel[2].line_of_sight.first_point.r == 6.2
        # Try appending new interfermeter channels
        add_interferometer!(
            "$(@__DIR__)/../samples/test_interferometer_new_channels.json", ids;
            n_e_gsi=-5,
        )
        for ch ∈ ids.interferometer.channel
            println()
            println("-"^49)
            println("-"^49)
            println("Channel $(ch.name)")
            println("-"^49)
            @printf("|%15s|%15s|%15s|\n", "time", "int_n_e", "n_e_average")
            println("-"^49)
            for ii ∈ 1:20:length(ch.n_e_line.time)
                @printf(
                    "|%15.3e|%15.3e|%15.3e|\n",
                    ch.n_e_line.time[ii],
                    ch.n_e_line.data[ii],
                    ch.n_e_line_average.data[ii]
                )
            end
            println("-"^49)
        end
        @test length(ids.interferometer.channel) == 5
        @test ids.interferometer.channel[5].name == "V1.5"
        @test ids.interferometer.channel[5].line_of_sight.first_point.r == 5.5
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
            n_e_gsi=-5,
        )
        # Just checking if the function runs through for now
        for lp ∈ ids.langmuir_probes.embedded
            println()
            println("-"^97)
            println("-"^97)
            println("Probe: $(lp.name)")
            println("-"^97)
            @printf("|%15s|%15s|%15s|%15s|%15s|%15s|\n", "time", "n_e", "t_e", "t_i",
                "i_sat", "j_sat")
            println("-"^97)
            for ii ∈ 1:20:length(lp.time)
                @printf(
                    "|%15.3e|%15.3e|%15.3e|%15.3e|%15.3e|%15.3e|\n",
                    lp.time[ii],
                    lp.n_e.data[ii],
                    lp.t_e.data[ii],
                    lp.t_i.data[ii],
                    lp.ion_saturation_current.data[ii],
                    lp.j_i_parallel.data[ii],
                )
            end
            println("-"^97)
        end
        @test true
        # Create an I-V curve and plot it
        v_probe = collect(-100:0.1:100)
        v_plasma_1 = 20.0
        v_plasma_2 = 0.0
        v_plasma_3 = -20.0
        Δv_1 = v_probe .- v_plasma_1
        Δv_2 = v_probe .- v_plasma_2
        Δv_3 = v_probe .- v_plasma_3
        Te = 10.0 # eV
        Ti = 9.0 # eV
        ne = 1e21 # m^-3
        A = 1e-6 # m^2 1 sq mm
        i_probe_1 = langmuir_probe_current.(Δv_1, Te, Ti, ne, A)
        i_probe_2 = langmuir_probe_current.(Δv_2, Te, Ti, ne, A)
        i_probe_3 = langmuir_probe_current.(Δv_3, Te, Ti, ne, A)
        plot(v_probe, i_probe_1; label="v_plasma = 20.0 V")
        plot!(v_probe, i_probe_2; label="v_plasma = 0.0 V")
        plot!(v_probe, i_probe_3; label="v_plasma = -20.0 V")
        plot!(; legend=true,
            xlabel="Probe Voltage (V)",
            ylabel="Probe Current (A)",
            title="Langmuir Probe I-V Curve",
        )
        details = text("Te = $Te eV\nTi = $Ti eV\nne = $ne m^-3\nA = $A m^2", :left)
        annotate!(-100, 300, details)
        savefig("$(@__DIR__)/langmuir_probe_iv.png")
        @test true
    end
end

function test_gas_response(config, excitation, plot_title, figname; fit=false)
    ids = IMAS.dd()
    add_gas_injection!(config, ids)
    ttotal = 5
    nt = Int(ttotal * 1000) + 1
    tstart = 1.0
    tend = 3.0
    toffset = 0.5
    tstartind = Int(tstart * 1000) + 1
    tendind = Int(tend * 1000) + 1
    tstartind2 = Int((tstart + toffset) * 1000) + 1
    tendind2 = Int((tend + toffset) * 1000) + 1

    tt = collect(LinRange(0, ttotal, nt))
    ids.gas_injection.valve[1].voltage.time = tt
    ids.gas_injection.valve[1].voltage.data = zeros(nt)
    ids.gas_injection.valve[1].voltage.data[tstartind:tendind] .=
        excitation.(tt[tstartind:tendind])

    ids.gas_injection.valve[2].voltage.time = tt
    ids.gas_injection.valve[2].voltage.data = zeros(nt)
    ids.gas_injection.valve[2].voltage.data[tstartind2:tendind2] .=
        excitation.(tt[tstartind:tendind])

    if fit
        # Setting special latency for GASD, GASA will follow global latency
        cmd_data = readdlm("$(@__DIR__)/../samples/gi193607_gacgasd.txt"; comments=true)
        P_ves_data =
            readdlm("$(@__DIR__)/../samples/gi193607_pcm240tor.txt"; comments=true)
        cmd_tt = cmd_data[:, 1] ./ 1000.0 # Convert to seconds
        cmd = cmd_data[:, 2]
        P_ves_tt = P_ves_data[:, 1] ./ 1000.0 # Convert to seconds
        P_ves = P_ves_data[:, 2] * 1e-3 * 133.3223684211 # Convert to Pa
        V_ves = 37.0 # m^3

        gasd_resp_curve, gasd_model =
            get_gas_injection_response(cmd, cmd_tt, P_ves, P_ves_tt, V_ves)
        ids.gas_injection.valve[1].response_curve = gasd_resp_curve
        # Adding made up time_constant and damping
        gasd_model[:time_constant] = 0.05
        gasd_model[:damping] = 0.5
        gasd_model[:dribble_decay_time_constant] = 0.45

        valves = Dict("GASD" => gasd_model)

        # Show a plot of how fitting worked.
        valve_response = linear_interpolation(
            gasd_resp_curve.voltage,
            gasd_resp_curve.flow_rate,
        )
        flow_rate = zeros(length(cmd_tt))
        tt0 = cmd_tt[1]
        tt_over_lat = findall(x -> x > gasd_model[:latency] + tt0, cmd_tt)
        if length(tt_over_lat) > 0
            skip = tt_over_lat[1]
            flow_rate[skip:end] = valve_response.(cmd[1:end-skip+1])
            flow_rate = map((x)::Float64 -> x < 0.0 ? 0.0 : x, flow_rate)
        end

        meas = P_ves * V_ves
        fitted = cumsum(flow_rate .* [0.0; diff(cmd_tt)])

        plot(P_ves_tt, meas; label="Measured")
        plot!(cmd_tt, fitted; label="Fitted gas model", linewidth=2, linestyle=:dash)
        plot!(;
            legend=true,
            legendposition=:topleft,
            xlabel="Time (s)",
            ylabel="Total accumulated gas / Pa m^3",
            title="Gas Injection Calibration",
        )
        savefig("$(@__DIR__)/gasd_calibration_fit.png")
    else
        # Setting special latency for GASD, GASA will follow global latency
        valves = Dict{String, Dict{Symbol, Any}}(
            "GASD" =>
                Dict(:latency => 0.183, :time_constant => 0.05, :damping => 0.5,
                    :dribble_decay_time_constant => 0.45),
        )
    end

    compute_gas_injection!(ids; valves=valves)

    plot(;
        title=plot_title,
        legend=true,
        xlabel="Time (s)",
        ylabel="Flow rate or Command Voltage",
    )
    for valve ∈ ids.gas_injection.valve
        plot!(
            valve.voltage.time,
            valve.voltage.data;
            lw=2,
            alpha=0.5,
            label="$(valve.name): Command Voltage / (Pa m^3/s)",
        )

        plot!(
            valve.flow_rate.time,
            valve.flow_rate.data;
            lw=2,
            label="$(valve.name): Flow Rate / (Pa m^3/s)",
        )
    end
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
            "$(@__DIR__)/gas_injection_step.png";
            fit=true,
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

if args["magic"]
    @testset "magic_diagnostics" begin
        ids = json2imas(
            "$(@__DIR__)/../samples/time_dep_edge_profiles_with_equilibrium.json",
        )
        nt = length(ids.edge_profiles.ggd)
        nesep = magic_nesep(ids; cell_grid_subset=-5)
        @test length(nesep) == nt
    end
end

if args["derived"]
    @testset "derived_quantities_simple" begin
        # Test data for power calculations
        time = [i for i ∈ 0:0.01:10.0]
        Wflat = 2.5e6  # J
        tau = 0.325  # s
        ramp_up_end = 1.15
        ramp_down_start = 9.0
        W = (time ./ ramp_up_end) .* Wflat
        W[time.>ramp_up_end] .= Wflat
        sel = time .> ramp_down_start
        W[sel] .= Wflat * (1 .- (time[sel] .- ramp_down_start))
        sel = (time .> 5.5) .& (time .< 6.8)
        W[sel] .+= Wflat .* 0.12 .* sin.(time[sel] .* 2 .* π ./ 0.5)
        dWdt = gradient(time, W)
        test_time = 0.9
        test_dWdt = (dWdt[abs.(time .- test_time).<0.011])[1]

        Pscale = Wflat / tau
        P_rad_core = Pscale * 0.32
        P_input = Pscale + P_rad_core
        P_NBI = P_input * 0.5
        P_ECH = P_input * 0.5
        P_OHM = P_input - (P_NBI + P_ECH)

        nt = length(time)

        P_SOL = calc_loss_power(test_dWdt, P_rad_core, P_OHM; P_NBI=P_NBI, P_ECH=P_ECH)
        @test P_SOL > 0

        tau_arr = 0.1 .+ time .* 0.0
        tau_arr[W.>1e6] .= tau
        Pscale_arr = W ./ tau_arr
        P_rad_core_arr = Pscale_arr .* 0.32
        P_input_arr = Pscale_arr .+ P_rad_core_arr
        P_NBI_arr = P_input_arr .* 0.5
        P_ECH_arr = P_input_arr .* 0.4
        P_OHM_arr = P_input_arr .- (P_NBI_arr .+ P_ECH_arr)

        P_SOL_arr = calc_loss_power(
            time,
            W,
            P_rad_core_arr,
            P_OHM_arr;
            P_NBI=P_NBI_arr,
            P_ECH=P_ECH_arr,
        )
        @test length(P_SOL_arr) == nt

        P_cond = calc_conducted_loss_power(test_dWdt, P_OHM; P_NBI=P_NBI, P_ECH=P_ECH)
        @test P_cond > 0

        P_cond_arr = calc_conducted_loss_power(
            time,
            W,
            P_OHM_arr;
            P_NBI=P_NBI_arr,
            P_ECH=P_ECH_arr,
        )
        @test length(P_cond_arr) == nt

        # Test data for qcyl
        B_ϕ_axis = -2.07  # T
        Iₚ = 1.1e6  # A
        aₘᵢₙₒᵣ = 0.56  # m
        R_geo = 1.675  # m
        κ = 1.81  # unitless
        q_cyl = calc_q_cyl(B_ϕ_axis, Iₚ, aₘᵢₙₒᵣ, R_geo, κ)
        println(q_cyl)
        @test q_cyl != 0

        λq1 = calc_heat_flux_width(B_ϕ_axis, P_SOL * 1e-6, R_geo, aₘᵢₙₒᵣ, κ, Iₚ)  # mm
        λq2 = calc_heat_flux_width(B_ϕ_axis, P_SOL * 1e-6, R_geo, q_cyl)  # mm
        @test λq1 == λq2

        B_θ_OMP = -0.2545  # T
        λq3 = calc_heat_flux_width(B_θ_OMP)  # mm
        @test λq3 > 0
    end
    @testset "Derived quantities from DD" begin
        ids = json2imas(
            "$(@__DIR__)/../samples/time_dep_edge_profiles_with_equilibrium.json",
        )
        R_OMP1, Z_OMP1 = find_OMP_RZ(ids)
        R_OMP, Z_OMP = find_OMP_RZ(ids.edge_profiles.grid_ggd[1])
        @test R_OMP > 0
        @test abs(Z_OMP) < R_OMP
        @test R_OMP1 == R_OMP
        @test Z_OMP1 == Z_OMP

        ids = json2imas("$(@__DIR__)/../samples/sd_input_data.json")
        B_θ_OMP = read_B_theta_OMP(
            ids;
            grid_ggds=ids.edge_profiles.grid_ggd,
            cell_grid_subset=-5,
        )
        @test length(B_θ_OMP) == length(ids.equilibrium.time_slice)

        nt = length(ids.equilibrium.time_slice)
        # ////////// TEMPORARY: this should go into the sample file later
        for eqt ∈ ids.equilibrium.time_slice
            p2 = eqt.profiles_2d[1]
            p2.grid_type.index = 1  # 1 = rectangular, such as dim1 = R, dim2 = Z
            p2.grid_type.name = "R-Z grid for flux map"
            p2.grid_type.description = (
                "A recntangular grid of points in R,Z on which poloidal " *
                "magnetic flux psi is defined. The grid's dim1 is R, dim2 is Z."
            )
        end
        ids.summary.global_quantities.power_ohm.value = zeros(nt) .+ 0.457
        ids.summary.heating_current_drive.power_nbi.value = zeros(nt) .+ 1.56782
        ids.summary.heating_current_drive.power_ec.value = zeros(nt) .+ 2.8
        ids.summary.heating_current_drive.power_ic.value = zeros(nt)
        ids.summary.heating_current_drive.power_lh.value = zeros(nt)
        ids.summary.heating_current_drive.power_additional.value = zeros(nt)
        ids.summary.boundary.geometric_axis_r.value = zeros(nt) .+ 6.285
        ids.summary.fusion.power.value = zeros(nt)
        # ///////////// END TEMPORARY data that is being added to sample until it can be updated
        IMAS.flux_surfaces(ids.equilibrium)  # provides volume and stuff
        summarize_flux_surfaces!(ids)
        ids.summary.global_quantities.power_radiated_inside_lcfs.value = zeros(nt)  # This is needed but won't be immediately available in the sample

        P_SOL = calc_loss_power(ids)
        nt = length(ids.summary.time)
        @test length(P_SOL) == nt

        λq1 = calc_heat_flux_width(ids; version=1)
        λq2 = calc_heat_flux_width(ids; version=2)
        @test length(λq1) == nt
        @test length(λq2) == nt

        B_θ_OMP_no_ggd = zeros(length(ids.equilibrium.time_slice))
        for time_idx ∈ 1:length(ids.equilibrium.time_slice)
            B_θ_OMP_no_ggd[time_idx] = read_B_theta_OMP_no_ggd(ids; time_idx=time_idx)
        end
        println(B_θ_OMP_no_ggd)
    end
end
