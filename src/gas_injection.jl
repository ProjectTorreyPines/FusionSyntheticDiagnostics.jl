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
            println("Computing flow rate for valve $(valve.name)")
            start_ind = length(valve.flow_rate.data) + 1
            end_ind = length(valve.voltage.data)
            t = valve.voltage.time
            cmd = valve.voltage.data
            if haskey(valves[vid], "dribbling_model")
                dbm = valves[vid]["dribbling_model"]
                println("Using dribbling_model")
                if length(valve.flow_rate.data) == 0
                    # Initialize dribbling model state variables
                    for key ∈ ["Gamma_main"]
                        dbm[key] = []
                    end
                    valve.flow_rate.data = [0.0]
                    valve.flow_rate.time = [t[1]]
                    start_ind = 2
                    dbm["t0"] = Inf
                    Gamma_end = (t) -> 0.0
                end

                if dbm["t0"] == Inf
                    for ii ∈ start_ind:end_ind
                        if cmd[ii] > 0
                            dbm["t0"] = valve.voltage.time[ii] # First time valve opens
                            break
                        end
                    end
                end
                Gamma_spike =
                    (t) ->
                        if t > dbm["t0"] + dbm["delta_on"]
                            dt = dbm["t0"] + dbm["delta_on"] - t
                            return dbm["Gamma_s"] * exp(dt / dbm["tau_spike"])
                        else
                            return 0.0  # Do nothing before command time + on delay
                        end
                for ii ∈ start_ind:end_ind
                    command_step = cmd[ii] - cmd[ii-1]
                    append!(
                        dbm["Gamma_main"],
                        [get_Gamma_main(t[ii], command_step, dbm)],
                    )
                    total_Gamma_main = cumulate_Gamma(dbm["Gamma_main"])
                    Gamma_end = update_Gamma_end(
                        Gamma_end,
                        t[ii],
                        cmd[ii-1:ii],
                        dbm)
                    append!(
                        valve.flow_rate.data,
                        Gamma_spike(t[ii]) + total_Gamma_main(t[ii]) +
                        Gamma_end(t[ii]),
                    )
                    append!(valve.flow_rate.time, t[ii])
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
                        valve_response(cmd[ii]),
                    )
                    append!(valve.flow_rate.time, t[ii])
                end
            end
        end
    end
end

function instant_gas_model(command, config)
    return config["p1"] .* (sqrt.(((command .* config["p2"]) .^ 2.0 .+ 1) .- 1))
end

function get_Gamma_main(command_time, command_step, dbm)
    Gamma_0 = instant_gas_model(command_step, dbm) * sign(command_step)
    function Gamma_main(t)
        if t > command_time + dbm["delta_on"] && Gamma_0 != 0.0
            dt = command_time + dbm["delta_on"] - t
            decay_on = (1 - exp(dt / dbm["tau_on"]))
            decay_over = (1 + dbm["Gamma_over_0"] * exp(dt / dbm["tau_over"]))
            return Gamma_0 * decay_on * decay_over
        else
            return 0.0 # Do nothing before command time + on delay
        end
    end
    return Gamma_main
end

function update_Gamma_end(
    Gamma_end,
    command_time,
    command_recent_two,
    dbm,
)
    if command_recent_two[1] > 0 && command_recent_two[2] == 0
        Gamma_main_array = copy(dbm["Gamma_main"])
        cum_Gamma_main = cumulate_Gamma(Gamma_main_array)
        Gamma_end =
            (t) -> begin
                start_main_val = cum_Gamma_main(command_time + dbm["delta_on"])
                if t < command_time + dbm["delta_on"]
                    return 0.0 # Do nothing before command time + on delay
                elseif t >= command_time + dbm["delta_on"] + dbm["D"] # End2 case
                    last_end1_val = start_main_val * exp(-dbm["D"] / dbm["tau_off1"])
                    dt = command_time + dbm["delta_on"] + dbm["D"] - t
                    return last_end1_val * exp(dt / dbm["tau_off2"]) - cum_Gamma_main(t)
                else # End1 case
                    dt = command_time + dbm["delta_on"] - t
                    return start_main_val * exp(dt / dbm["tau_off1"]) -
                           cum_Gamma_main(t)
                end
            end
    elseif command_recent_two[2] > 0 && command_recent_two[1] == 0
        # Return identical zero if valve is turned back on
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