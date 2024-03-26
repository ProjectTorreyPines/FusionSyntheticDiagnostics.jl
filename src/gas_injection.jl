import Interpolations: linear_interpolation
import DSP.Filters: Biquad, bilinear, convert, SecondOrderSections, DF2TFilter
import DSP: filt

default_gas_injection = "$(@__DIR__)/default_gas_injection.json"

"""
    add_gas_injection!(
    config::String=default_gas_injection,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite::Bool=false,
    verbose::Bool=false,

)::IMASDD.dd

Add gas valves from a JSON file and compute the gas flow rate based on the command
signal in the gas valves.
"""
function add_gas_injection!(
    config::String=default_gas_injection,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite::Bool=false,
    verbose::Bool=false,
)::IMASDD.dd
    if endswith(config, ".json")
        config_dict = convert_strings_to_symbols(IMASDD.JSON.parsefile(config))
        add_gas_injection!(config_dict, ids; overwrite=overwrite, verbose=verbose)
    else
        error("Only JSON files are supported.")
    end
    return ids
end

"""
    add_gas_injection!(
    config::Dict{Symbol, Any},
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite::Bool=false,
    verbose::Bool=false,

)::IMASDD.dd

Add gas valves from a dictionary and compute the gas flow rate based on the command
signal in the gas valves.
"""
function add_gas_injection!(
    config::Dict{Symbol, Any},
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite::Bool=false,
    verbose::Bool=false,
)::IMASDD.dd
    # Check for duplicates
    if length(ids.gas_injection.valve) > 0
        duplicate_indices = []
        new_valves = [ch[:name] for ch ∈ config[:gas_injection][:valve]]

        for (ii, ch) ∈ enumerate(ids.gas_injection.valve)
            if ch.name ∈ new_valves
                append!(duplicate_indices, ii)
            end
        end
        if overwrite
            for ii ∈ reverse(duplicate_indices)
                println(
                    "Overwriting gas_injection valve ",
                    "$(ids.gas_injection.valve[ii].name)...",
                )
                deleteat!(ids.gas_injection.valve, ii)
            end
        else
            if length(duplicate_indices) > 0
                err_msg =
                    "Duplicate gas_injection valves found with " *
                    "overlapping names.\n"
                for ii ∈ duplicate_indices
                    err_msg *= "$(ids.gas_injection.valve[ii].name)\n"
                end
                err_msg *= "Use overwrite=true to replace them."
                throw(OverwriteAttemptError(err_msg))
            end
        end
        config[:gas_injection] =
            mergewith(
                append!,
                IMASDD.imas2dict(ids.gas_injection),
                config[:gas_injection],
            )
    end
    IMASDD.dict2imas(config, ids; verbose=verbose)
    valves = Dict{String, Dict{Symbol, Any}}(
        valve[:name] => valve for valve ∈ config[:gas_injection][:valve]
    )
    compute_gas_injection!(ids; valves=valves)
    return ids
end

"""
    compute_gas_injection!(ids::IMASDD.dd)

Compute the gas flow rate based on the command signal in the gas valves.
"""
function compute_gas_injection!(
    ids::IMASDD.dd;
    valves::Dict{String, Dict{Symbol, Any}}=Dict{String, Dict{Symbol, Any}}(),
)
    if IMASDD.ismissing(ids.gas_injection, :latency)
        global_latency = 0.0
    elseif IMASDD.isempty(ids.gas_injection.latency)
        global_latency = 0.0
    else
        global_latency = ids.gas_injection.latency
    end

    for valve ∈ ids.gas_injection.valve
        proceed =
            !IMASDD.ismissing(valve.response_curve, :flow_rate) &&
            !IMASDD.ismissing(valve.response_curve, :voltage) &&
            !IMASDD.ismissing(valve.voltage, :data) &&
            !IMASDD.ismissing(valve.voltage, :time)
        if proceed
            if !IMASDD.isempty(valve.response_curve.flow_rate) &&
               !IMASDD.isempty(valve.response_curve.voltage)
                if length(valve.response_curve.flow_rate) !=
                   length(valve.response_curve.voltage)
                    error(
                        "$(valve.name): Length of flow rate and voltage in " *
                        "response_curve should be equal.",
                    )
                end
            else
                proceed = false
                println(
                    "Response curve for valve $(valve.name):  is empty. " \
                    "Skipping computation of flow rate.",
                )
            end
            if !IMASDD.isempty(valve.voltage.data) &&
               !IMASDD.isempty(valve.voltage.time)
                if length(valve.voltage.data) != length(valve.voltage.time)
                    error(
                        "$(valve.name): Length of data and time in voltage " *
                        "should be equal.",
                    )
                end
            else
                proceed = false
                println(
                    "Voltage data for valve $(valve.name) is empty. " \
                    "Skipping computation of flow rate.",
                )
            end
        end
        if proceed
            latency = deepcopy(global_latency)
            LPF = nothing
            if valve.name ∈ keys(valves)
                valve_model = valves[valve.name]
                if :latency ∈ keys(valve_model)
                    latency = valve_model[:latency]
                end
                if :time_constant ∈ keys(valve_model) && :damping ∈ keys(valve_model)
                    LPF = get_lpf(
                        1 / (valve.voltage.time[2] - valve.voltage.time[1]),
                        valve_model[:time_constant],
                        valve_model[:damping],
                        1.0,
                    )
                end
            end
            valve.flow_rate.time = valve.voltage.time
            flow_rate = zeros(length(valve.voltage.time))
            valve_response = linear_interpolation(
                valve.response_curve.voltage,
                valve.response_curve.flow_rate,
            )
            tt0 = valve.voltage.time[1]
            tt_over_lat = findall(x -> x > latency + tt0, valve.voltage.time)
            if length(tt_over_lat) > 0
                skip = tt_over_lat[1]
                flow_rate[skip:end] = valve_response.(valve.voltage.data[1:end-skip+1])
                if !isnothing(LPF)
                    flow_rate = filt(LPF, flow_rate)
                end
                flow_rate = map((x)::Float64 -> x < 0.0 ? 0.0 : x, flow_rate)
            end
            valve.flow_rate.data = flow_rate
        end
    end
end

"""
    get_lpf(fs::Float64, tau::Float64, damping::Float64, gain::Float64)

Create a second order filter with the given parameters. The filter is created in the
s-domain and then converted to the z-domain using bilinear transform.
In s-domain, the filter transfer function is:

    H(s) = gain * ωₙ^2 / (s^2 + 2 * damping * ωₙ * s + ωₙ^2)

where ωₙ = 2π / tau
"""
function get_lpf(fs::Float64, tau::Float64, damping::Float64, gain::Float64)
    ωₙ = 2π / tau
    b = [0.0, 0.0, ωₙ^2] .* gain
    a = [2 * damping * ωₙ, ωₙ^2]
    filter_s = Biquad{:s}(b..., a...)
    return convert(SecondOrderSections, bilinear(filter_s, fs))
end