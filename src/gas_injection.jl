using JSON: JSON
import Interpolations: linear_interpolation

default_gas_injection = "$(@__DIR__)/default_gas_injection.json"

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

function compute_gas_injection!(ids::OMAS.dd, valves::Vector{Any})
    for vid ∈ eachindex(ids.gas_injection.valve)
        valve = ids.gas_injection.valve[vid]
        if length(valve.voltage.data) > length(valve.flow_rate.data)
            start_ind = length(valve.flow_rate) + 1
            end_ind = length(valve.voltage.command)
            if haskey(valves[vid], "dribbling_model")
                dbm = valves[vid]["dribbling_model"]
                if length(valve.flow_rate.data) == 0
                    # Initialize dribbling model state variables
                    for key ∈ ["Gamma_spike", "Gamma_main", "Gamma_end"]
                        dbm[key] = []
                    end
                    # for key ∈ ["t0", "tf"]
                    #     dbm[key] = Inf
                    # end
                end

                # if dbm["t0"] == Inf
                #     for ii ∈ start_ind:end_ind
                #         if valve.voltage.command[ii] > 0
                #             dbm["t0"] = valve.voltage.time[ii] # First time valve opens
                #             break
                #         end
                #     end
                # end
                # if dbm["tf"] == Inf && dbm["t0"] != Inf
                #     for ii ∈ start_ind:end_ind
                #         if valve.voltage.command[ii] == 0 && valve.voltage.time[ii] > dbm["t0"]
                #             dbm["tf"] = valve.voltage.time[ii] # First time valve closes
                #             break
                #         end
                #     end
                # end
                for ii ∈ start_ind:end_ind
                    append!(
                        dbm["Gamma_spike"],
                        get_Gamma_spike(t[ii], valve.voltage.command[ii-1:ii], dbm),
                    )
                    command_step = valve.voltage.data[ii] - valve.voltage.data[ii-1]
                    append!(dbm["Gamma_main"], get_Gamma_main(t[ii], command_step, dbm))
                    append!(
                        dbm["Gamma_end"],
                        get_Gamma_end(
                            t[ii],
                            valve.voltage.command[ii-1:ii],
                            dbm,
                            cumulate_Gamma(dbm["Gamma_main"]),
                        ),
                    )
                    total_Gamma_spike = cumulate_Gamma(dbm["Gamma_spike"])
                    total_Gamma_main = cumulate_Gamma(dbm["Gamma_main"])
                    total_Gamma_end = cumulate_Gamma(dbm["Gamma_end"])
                    append!(
                        valve.flow_rate.data,
                        total_Gamma_spike(t[ii]) + total_Gamma_main(t[ii]) +
                        total_Gamma_end(t[ii]),
                    )
                end
            elseif length(valve.response_curve.flow_rate) ==
                   length(valve.response_curve.voltage) > 1
                valve_response = linear_interpolation(
                    valve.response_curve.voltage,
                    valve.response_curve.flow_rate,
                )
                for ii ∈ start_ind:end_ind
                    append!(
                        valve.flow_rate.data,
                        valve_response(valve.voltage.data[ii]),
                    )
                end
            end
        end
    end
end

function get_Gamma_spike(command_time, command_recent_two, dbm)
    if command_recent_two[1] == 0 && command_recent_two[2] > 0
        Gamma_spike = (t) ->
            if t > command_time + dbm["delta_on"]
                dt = command_time + dbm["delta_on"] - t
                return dbm["Gamma_s"] * exp(dt / dbm["tau_spike"])
            else
                return 0.0  # Do nothing before command time + on delay
            end
    else # Return identical zero if no turn on detected
        Gamma_spike = (t) -> 0.0
    end
    return Gamma_spike
end

function instant_gas_model(command, config)
    return config["p1"] .* (sqrt.(((command * config["p2"]) .^ 2.0 .+ 1) .- 1))
end

function get_Gamma_main(command_time, command_step, dbm)
    Gamma_0 = instant_gas_model(command_step, dbm)
    function Gamma_main(t)
        if t > command_time + dbm["delta_on"]
            dt = command_time + dbm["delta_on"] - t
            decay_on = (1 - exp(dt / dbm["tau_on"]))
            decay_over = (1 + dbm["Gamma_over_0"] / Gamma_0 * exp(dt / dbm["tau_over"]))
            return Gamma_0 * decay_on * decay_over
        else
            return 0.0 # Do nothing before command time + on delay
        end
    end
    return Gamma_main
end

function get_Gamma_end(command_time, command_recent_two, dbm, cum_Gamma_main)
    if command_recent_two[1] > 0 && command_recent_two[2] == 0
        Gamma_end =
            (t) -> begin
                start_main_val = cum_Gamma_main(command_time + dbm["delta_on"])
                if t < command_time + dbm["delta_on"]
                    return 0.0 # Do nothing before command time + on delay
                elseif t >= command_time + dbm["delta_on"] + dbm["D"] # End1 case
                    last_end1_val = start_main_val * exp(-dbm["D"] / dbm["tau_spike"])
                    dt = command_time + dbm["delta_on"] + dbm["D"] - t
                    return last_end1_val * exp(dt / dbm["tau_off2"]) - cum_Gamma_main(t)
                else # End2 case
                    dt = command_time + dbm["delta_on"] - t
                    return start_main_val * exp(dt / dbm["tau_off1"]) -
                           cum_Gamma_main(t)
                end
            end
    else # Return identical zero if no turn off detected
        Gamma_end = (t) -> 0.0
    end
    return Gamma_end
end

function cumulate_Gamma(Gamma_array)
    function cum_Gamma(t)
        sum = 0.0
        for Gamma ∈ Gamma_array
            sum += Gamma(t)
        end
        return sum
    end
    return cum_Gamma
end