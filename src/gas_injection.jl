using JSON: JSON
import Interpolations: linear_interpolation
import DSP.Filters: Biquad, bilinear, convert, SecondOrderSections, DF2TFilter
import DSP: filt

default_gas_injection = "$(@__DIR__)/default_gas_injection.json"

"""
    add_gas_injection!(
    config::String=default_gas_injection,
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),

)::OMAS.dd

Add gas valves from a JSON file and compute the gas flow rate based on the command
signal in the gas valves.
"""
function add_gas_injection!(
    config::String=default_gas_injection,
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
)::OMAS.dd
    if endswith(config, ".json")
        OMAS.json2imas(config, ids)
    else
        error("Only JSON files are supported.")
    end
    json_data = JSON.parsefile(config)
    valves = json_data["gas_injection"]["valve"]
    compute_gas_injection!(ids, valves)
    return ids
end

"""
    compute_gas_injection!(ids::OMAS.dd, valves::Vector{Any})

Compute gas flow rate based on command signal in the gas valves.
Currently, it takes an additional argument `valves` which is a vector of dictionaries
containing the valve parameters read from the JSON file where the complete OMAS
structure is present. This is because the gas injection dribbling model is not part of
the OMAS datastructure yet. Once it is added, this function will be modified to take
only the OMAS datastructure as an argument and get valves from inside it.
"""
function compute_gas_injection!(ids::OMAS.dd, valves::Vector{Any})
    for vid ∈ eachindex(ids.gas_injection.valve)
        valve = ids.gas_injection.valve[vid]
        if length(valve.voltage.data) > length(valve.flow_rate.data)
            println("Computing flow rate for valve $(valve.name)")
            t = valve.voltage.time
            cmd = valve.voltage.data
            if haskey(valves[vid], "dribbling_model")
                dbm = valves[vid]["dribbling_model"]
                println("Using dribbling_model")

                # Find time segments of on and off stages
                time_seg_inds = []
                last_start = 1
                for ii ∈ 2:length(cmd)
                    if cmd[ii] == 0 && cmd[ii-1] > 0 || cmd[ii] > 0 && cmd[ii-1] == 0
                        append!(time_seg_inds, [(last_start, ii - 1)])
                        last_start = ii
                    end
                end
                append!(time_seg_inds, [(last_start, length(cmd))])

                dead_time_shift = round(Int64, dbm["delta_on"] / (t[2] - t[1]))
                valve.flow_rate.data = zeros(Float64, dead_time_shift)
                popped = [false]
                for ts ∈ time_seg_inds
                    t_seg = t[ts[1]:ts[2]]
                    c_seg = cmd[ts[1]:ts[2]]
                    if ts[1] == 1
                        start_val = 0.0
                    else
                        start_val = valve.flow_rate.data[end]
                    end
                    if cmd[ts[1]] == 0
                        append!(
                            valve.flow_rate.data,
                            get_Gamma_end(dbm, t_seg, start_val),
                        )
                    else
                        append!(
                            valve.flow_rate.data,
                            get_Gamma_main(dbm, t_seg, c_seg, start_val, popped),
                        )
                    end
                end
                valve.flow_rate.data = valve.flow_rate.data[1:length(cmd)]
                valve.flow_rate.time = t
            elseif length(valve.response_curve.flow_rate) ==
                   length(valve.response_curve.voltage) > 1
                valve_response = linear_interpolation(
                    valve.response_curve.voltage,
                    valve.response_curve.flow_rate,
                )
                for ii ∈ eachindex(cmd)
                    append!(
                        valve.flow_rate.data,
                        valve_response(cmd[ii]),
                    )
                    append!(valve.flow_rate.time, t[ii])
                end
            end
        end
    end
end

"""
    instant_gas_model(command, config)

Ideal gas flow rate model that instantaneously changes the flow rate based on the
command voltage. The model is based on the following equation:

    Gamma = p1 * sqrt((p2 * command)^2 + 1) - p1

where p1 (Pa m³/s)and p2 (V⁻¹) are parameters of the model.
"""
function instant_gas_model(command, config)
    return config["p1"] .* (sqrt.(((command .* config["p2"]) .^ 2.0 .+ 1) .- 1))
end

"""
    get_Gamma_end(dbm, time_seg, start_val)

Get the flow rate for the end stage of the valve operation. This is the part where
the valve is closed and the gas is still flowing out of the valve until the gas
pressure drops to very low value. This happens in two stages:

 1. The gas pressure drops exponentially but very slowly, and it looks almost linear.
    That is because this stage only occurs for a short time "D" taken from dbm model
    dictionary.
 2. The gas pressure drops very fast. This is the part where the gas pressure drops to
    zero and flow of gas stops.
"""
function get_Gamma_end(dbm, time_seg, start_val)
    if time_seg[end] - time_seg[1] < dbm["D"]
        return start_val * exp.((time_seg[1] .- time_seg) / dbm["tau_off1"])
    else
        change_ind = findfirst(time_seg .> time_seg[1] + dbm["D"])
        t1 = time_seg[1:change_ind-1]
        t2 = time_seg[change_ind:end]
        ret = start_val * exp.((time_seg[1] .- t1) / dbm["tau_off1"])
        ret = append!(ret, ret[end] * exp.((t1[end] .- t2) / dbm["tau_off2"]))
        return ret
    end
end

"""
    get_filter(fs, tau, damping, gain)

Create a second order filter with the given parameters. The filter is created in the
s-domain and then converted to the z-domain using bilinear transform.
In s-domain, the filter transfer function is:

    H(s) = gain * ωₙ^2 / (s^2 + 2 * damping * ωₙ * s + ωₙ^2)

where ωₙ = 2π / tau
"""
function get_filter(fs, tau, damping, gain)
    ωₙ = 2π / tau
    b = [0.0, 0.0, ωₙ^2] .* gain
    a = [2 * damping * ωₙ, ωₙ^2]
    filter_s = DSP.Filters.Biquad{:s}(b..., a...)
    return convert(SecondOrderSections, bilinear(filter_s, fs))
    # return DF2TFilter(filter_z, si)
end

"""
    get_Gamma_main(dbm, t_seg, c_seg, start_val, popped)

Get the flow rate for the main stage of the valve operation. If popped is [false], it
will add a spike at the begining of flow rate that happens due to valve being
popped open with a large short pulse.
Note that popped is a boolean array of size 1. This is because Julia passes arrays by
reference, so we can modify the value of popped in the function. Once popped, it won't
do it again.
"""
function get_Gamma_main(dbm, t_seg, c_seg, start_val, popped)
    fs = 1 / (t_seg[2] - t_seg[1])
    main_filter = get_filter(
        fs,
        dbm["tau_on"],
        dbm["damping_on"],     # oscillating < 0.707 < over-damped
        dbm["p1"] * dbm["p2"], # This takes the linear part of ideal gas injection
    ) # Returns SOS filter for main valve operation, no initial state
    si = zeros(Float64, 2, 1)
    si[1] = start_val / main_filter.g - main_filter.biquads[1].b0 * c_seg[1]
    main_filter_sf = DF2TFilter(main_filter, si) # Stateful filter with start value
    if popped[1]
        return filt(main_filter_sf, c_seg)
    else
        spike = dbm["Gamma_s"] * exp.((t_seg[1] .- t_seg) ./ dbm["tau_spike"])
        popped[1] = true
        return filt(main_filter_sf, c_seg) + spike
    end
end
