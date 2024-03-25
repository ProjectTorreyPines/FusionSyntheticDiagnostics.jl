import Interpolations: linear_interpolation
# import DSP.Filters: Biquad, bilinear, convert, SecondOrderSections, DF2TFilter
# import DSP: filt

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
    compute_gas_injection!(ids)
    return ids
end

"""
    compute_gas_injection!(ids::IMASDD.dd)

Compute the gas flow rate based on the command signal in the gas valves.
"""
function compute_gas_injection!(ids::IMASDD.dd)
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
            valve.flow_rate.time = valve.voltage.time
            valve_response = linear_interpolation(
                valve.response_curve.voltage,
                valve.response_curve.flow_rate,
            )
            valve.flow_rate.data = valve_response.(valve.voltage.data)
        end
    end
end
