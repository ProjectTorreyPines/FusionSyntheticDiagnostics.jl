import Interpolations: linear_interpolation, Flat
import DSP.Filters: Biquad, bilinear, convert, SecondOrderSections, DF2TFilter
import DSP: filt, xcorr
import LsqFit: curve_fit, coef
import Statistics: mean

export add_gas_injection!, compute_gas_injection!, compute_gas_injection,
    get_gas_injection_response, get_required_gas_cmd

default_gas_injection = "$(@__DIR__)/default_gas_injection.json"

"""
    add_gas_injection!(
        config::String=default_gas_injection,
        @nospecialize(ids::IMAS.dd)=IMAS.dd();
        overwrite::Bool=false,
        verbose::Bool=false,
    )::IMAS.dd

Add gas valves from a JSON file and compute the gas flow rate based on the command
signal in the gas valves.
"""
function add_gas_injection!(
    config::String=default_gas_injection,
    @nospecialize(ids::IMAS.dd)=IMAS.dd();
    overwrite::Bool=false,
    verbose::Bool=false,
)::IMAS.dd
    if endswith(config, ".json")
        config_dict = convert_strings_to_symbols(IMAS.JSON.parsefile(config)) # Use with import IMASdd as IMAS
        # config_dict = convert_strings_to_symbols(IMAS.IMASdd.JSON.parsefile(config)) # Use with using IMAS: IMAS
        add_gas_injection!(config_dict, ids; overwrite=overwrite, verbose=verbose)
    else
        error("Only JSON files are supported.")
    end
    return ids
end

"""
    add_gas_injection!(
        config::Dict{Symbol, Any},
        @nospecialize(ids::IMAS.dd)=IMAS.dd();
        overwrite::Bool=false,
        verbose::Bool=false,
    )::IMAS.dd

Add gas valves from a dictionary and compute the gas flow rate based on the command
signal in the gas valves.
"""
function add_gas_injection!(
    config::Dict{Symbol, Any},
    @nospecialize(ids::IMAS.dd)=IMAS.dd();
    overwrite::Bool=false,
    verbose::Bool=false,
)::IMAS.dd
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
                IMAS.imas2dict(ids.gas_injection),
                config[:gas_injection],
            )
    end
    IMAS.dict2imas(config, ids; verbose=verbose)
    valves = Dict{String, Dict{Symbol, Any}}(
        valve[:name] => valve for valve ∈ config[:gas_injection][:valve]
    )
    compute_gas_injection!(ids; valves=valves)
    return ids
end

"""
    compute_gas_injection!(
        ids::IMAS.dd;
        valves::Dict{String, Dict{Symbol, Any}}=Dict{String, Dict{Symbol, Any}}(),
    )::Array{Vector{Float64}}

Compute the gas flow rate based on the command signal in the all gas valves.
"""
function compute_gas_injection!(
    ids::IMAS.dd{T};
    valves::Dict{String, Dict{Symbol, Any}}=Dict{String, Dict{Symbol, Any}}(),
)::Array{Vector{Float64}} where {T}
    if IMAS.ismissing(ids.gas_injection, :latency)
        global_latency = 0.0
    elseif IMAS.isempty(ids.gas_injection.latency)
        global_latency = 0.0
    else
        global_latency = ids.gas_injection.latency
    end

    future_flow_rates = Array{Vector{Float64}}(undef, length(ids.gas_injection.valve))
    for (vind, valve) ∈ enumerate(ids.gas_injection.valve)
        future_flow_rates[vind] = Float64[]
        proceed =
            !IMAS.ismissing(valve.response_curve, :flow_rate) &&
            !IMAS.ismissing(valve.response_curve, :voltage) &&
            !IMAS.ismissing(valve.voltage, :data) &&
            !IMAS.ismissing(valve.voltage, :time)
        if proceed
            if !IMAS.isempty(valve.response_curve.flow_rate) &&
               !IMAS.isempty(valve.response_curve.voltage)
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
            if !IMAS.isempty(valve.voltage.data) &&
               !IMAS.isempty(valve.voltage.time)
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
            if valve.name ∈ keys(valves)
                valve_model = valves[valve.name]
            else
                valve_model = Dict{Symbol, Any}()
            end
            future_flow_rates[vind] = compute_gas_injection!(
                valve;
                valve_model=valve_model,
                global_latency=global_latency,
            )
        end
    end
    return future_flow_rates
end

"""
    compute_gas_injection(
        tt::Vector{Float64},
        cmd_voltage::Vector{Float64},
        response_curve_voltage::Vector{Float64},
        response_curve_flow_rate::Vector{Float64};
        valve_model::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        global_latency::Float64=0.0,
    )::Tuple{Vector{Float64}, Vector{Float64}, Union{Nothing, Vector{Float64}}}

Lowest level function to compute gas flow rate based on the command voltage data and
response cruve data all provided in base Julia types.
"""
function compute_gas_injection(
    tt::Vector{Float64},
    cmd_voltage::Vector{Float64},
    response_curve_voltage::Vector{Float64},
    response_curve_flow_rate::Vector{Float64};
    valve_model::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    global_latency::Float64=0.0,
)::Tuple{Vector{Float64}, Vector{Float64}, Union{Nothing, Vector{Float64}}}
    latency = deepcopy(global_latency)
    LPF = nothing
    dribble_tau = nothing
    if :latency ∈ keys(valve_model)
        latency = valve_model[:latency]
    end
    if :time_constant ∈ keys(valve_model) && :damping ∈ keys(valve_model) &&
       length(tt) > 1
        LPF = get_lpf(
            1 / (tt[2] - tt[1]),
            valve_model[:time_constant],
            valve_model[:damping],
            1.0,
        )
    end
    if :dribble_decay_time_constant ∈ keys(valve_model)
        dribble_tau = valve_model[:dribble_decay_time_constant]
    end
    flow_rate = zeros(length(tt))
    valve_response = linear_interpolation(
        response_curve_voltage,
        response_curve_flow_rate;
        extrapolation_bc=Flat(),
    )
    tt_over_lat = findall(x -> x >= latency + tt[1], tt)
    if length(tt_over_lat) > 0
        skip = tt_over_lat[1]
        flow_rate = valve_response.(cmd_voltage)
        if !isnothing(dribble_tau) && length(tt) > 1
            flow_rate = dribble(
                flow_rate,
                dribble_tau,
                1 / (tt[2] - tt[1]),
            )
        end
        if !isnothing(LPF)
            flow_rate = filt(LPF, flow_rate, create_si(LPF, flow_rate[1]))
        end
        flow_rate = map((x)::Float64 -> x < 0.0 ? 0.0 : x, flow_rate)
        future_flow_rates = flow_rate[end-skip+2:end]
        flow_rate[1:skip-1] .= 0.0
        flow_rate[skip:end] = flow_rate[1:end-skip+1]
    else
        future_flow_rates = nothing
    end
    return tt, flow_rate, future_flow_rates
end

"""
    compute_gas_injection(
        tt::Vector{Float64},
        cmd_voltage::Vector{Float64},
        response_curve::IMAS.gas_injection__valve___response_curve;
        valve_model::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        global_latency::Float64=0.0,
    )::Tuple{Vector{Float64}, Vector{Float64}, Union{Nothing, Vector{Float64}}}

Convinience function format where response curve is provided as the ids type but
everything else is provided in base Julia types.
"""
function compute_gas_injection(
    tt::Vector{Float64},
    cmd_voltage::Vector{Float64},
    response_curve::IMAS.gas_injection__valve___response_curve;
    valve_model::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    global_latency::Float64=0.0,
)::Tuple{Vector{Float64}, Vector{Float64}, Union{Nothing, Vector{Float64}}}
    return compute_gas_injection(
        tt,
        cmd_voltage,
        response_curve.voltage,
        response_curve.flow_rate;
        valve_model=valve_model,
        global_latency=global_latency,
    )
end

"""
    compute_gas_injection(
        valve::IMAS.gas_injection__valve;
        valve_model::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        global_latency::Float64=0.0,
    )::Tuple{Vector{Float64}, Vector{Float64}, Union{Nothing, Vector{Float64}}}

Top most level function to compute gas flow rate for a single valve.
"""
function compute_gas_injection(
    valve::IMAS.gas_injection__valve;
    valve_model::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    global_latency::Float64=0.0,
)::Tuple{Vector{Float64}, Vector{Float64}, Union{Nothing, Vector{Float64}}}
    return compute_gas_injection(
        valve.voltage.time,
        valve.voltage.data,
        valve.response_curve;
        valve_model=valve_model,
        global_latency=global_latency,
    )
end

"""
    compute_gas_injection!(
        valve::IMAS.gas_injection__valve;
        valve_model::Dict{Symbol, Any}=Dict{Symbol, Any}(),
        global_latency::Float64=0.0,
    )::Union{Nothing, Vector{Float64}}

In-place version of [`compute_gas_injection`](@ref) function for a single valve.
"""
function compute_gas_injection!(
    valve::IMAS.gas_injection__valve;
    valve_model::Dict{Symbol, Any}=Dict{Symbol, Any}(),
    global_latency::Float64=0.0,
)::Union{Nothing, Vector{Float64}}
    tt, flow_rate, future_flow_rates = compute_gas_injection(
        valve;
        valve_model=valve_model,
        global_latency=global_latency,
    )
    valve.flow_rate.time = tt
    valve.flow_rate.data = flow_rate
    return future_flow_rates
end

"""
    get_lpf(
        fs::Float64,
        tau::Float64,
        damping::Float64,
        gain::Float64,
    )::SecondOrderSections

Create a second order filter with the given parameters. The filter is created in the
s-domain and then converted to the z-domain using bilinear transform.
In s-domain, the filter transfer function is:

``H(s) = k \\frac{ωₙ^2}{ s^2 + 2 \\gamma ωₙ s + ωₙ^2}``

where ωₙ = 2π / `tau`, ``\\gamma`` is `damping` and `k` is `gain`.
"""
function get_lpf(
    fs::Float64,
    tau::Float64,
    damping::Float64,
    gain::Float64,
)::SecondOrderSections
    ωₙ = 2π / tau
    b = [0.0, 0.0, ωₙ^2] .* gain
    a = [2 * damping * ωₙ, ωₙ^2]
    filter_s = Biquad{:s}(b..., a...)
    return convert(SecondOrderSections, bilinear(filter_s, fs))
end

"""
    create_si(filter, last_val)

Create initial state for a filter based on a provided last value. This assumes that
upto this point, there was a constant command on (called prev_xi in this code) and the
filter had reached equilibrium to a constant output which is the last_val provided.
"""
function create_si(filter::SecondOrderSections, last_val::Float64)
    n = length(filter.biquads)
    g = filter.g
    si = zeros(Float64, 2, n)
    for fi ∈ n:-1:1
        b0 = filter.biquads[fi].b0
        b1 = filter.biquads[fi].b1
        b2 = filter.biquads[fi].b2
        a1 = filter.biquads[fi].a1
        a2 = filter.biquads[fi].a2
        if fi == n
            prev_yi = last_val / g
        else
            prev_yi = last_val
        end
        prev_xi = prev_yi * (1 + a1 + a2) / (b0 + b1 + b2)
        si[2, fi] = b2 * prev_xi - a2 * prev_yi
        si[1, fi] = si[2, fi] + b1 * prev_xi - a1 * prev_yi
        last_val = prev_xi
    end
    return si
end

"""
    dribble(
        data::Vector{Float64},
        decay_time_constant::Float64,
        fs::Float64,
    )::Vector{Float64}

Function to model dribble effect when gas command falls too sharply due to remaining
gas in the pipe. The function compares change in `data`(`ΔD`) at every time step with

``ΔD < -\\frac{1}{τ f_s} e^{ - \\left( \\frac{1}{τ f_s} \\right) }``

where `τ`(`decay_time_constant`) is the decay time constant and ``f_s``(`f_s`) is the
sampling frequency. If the above condition is found true, the function forces `data` to
decay exponentially with the decay time constant untill the condition is false. The
function returns the modified `data`.
"""
function dribble(
    data::Vector{Float64},
    decay_time_constant::Float64,
    fs::Float64,
)::Vector{Float64}
    data_der = [diff(data); 0.0]
    ii = 1
    while ii < length(data_der) - 1
        if data_der[ii] <
           -exp(-1.0 / (decay_time_constant * fs)) / (decay_time_constant * fs)
            data[ii+1] = data[ii] * exp(-1.0 / (decay_time_constant * fs))
            data_der[ii] = data[ii+1] - data[ii]
            if ii < length(data_der) - 2
                data_der[ii+1] = data[ii+2] - data[ii+1]
            end
        end
        ii += 1
    end
    return data
end

"""
    downsample_smooth(
        tt::Vector{Float64},
        data::Vector{Float64},
        dt_res::Float64;
        default_start=0.0,
    )::Vector{Float64}

Downsample and smooth the data to a new time resolution. The time axis does not need to
be equally spaced in time. The function creates a new time axis with spacing given by
`dt_res` and then averages the data points that fall within the new time bin. If no data
points fall within the new time bin, the function uses the last data point to fill the
new time bin. `default_start` is used to fill first time bin if no data is present in
the first time bin.
"""
function downsample_smooth(
    tt::Vector{Float64},
    data::Vector{Float64},
    dt_res::Float64;
    default_start=0.0,
)::Vector{Float64}
    tt_res = collect(tt[1]:dt_res:tt[end])
    return downsample_smooth(tt, data, tt_res; default_start=default_start)
end

"""
    downsample_smooth(
        tt::Vector{Float64},
        data::Vector{Float64},
        tt_res::Vector{Float64};
        default_start=0.0,
    )::Vector{Float64}

Downsample and smooth the data to a new time resolution. The time axis does not need to
be equally spaced in time. The function uses `tt_res` as the new resampled time axis and
averages the data points that fall within the new time bin. If no data points fall
within the new time bin, the function uses the last data point to fill the new time bin.
`default_start` is used to fill first time bin if no data is present in the first time
bin. Note that `tt_res` need not be equally spaced in time either.
"""
function downsample_smooth(
    tt::Vector{Float64},
    data::Vector{Float64},
    tt_res::Vector{Float64};
    default_start=0.0,
)::Vector{Float64}
    data_res = zeros(length(tt_res))
    last_time = -Inf
    for ii ∈ eachindex(data_res)
        snap = data[(tt.>last_time).&(tt.<=tt_res[ii])]

        if length(snap) == 0
            if ii > 1
                data_res[ii] = data_res[ii-1]
            else
                data_res[ii] = default_start
            end
        else
            data_res[ii] = mean(snap)
        end
        last_time = tt_res[ii]
    end
    return data_res
end

"""
    find_delay(
        data1::Vector{Float64},
        data2::Vector{Float64},
        dt::Float64,
    )::Float64

Find the delay between two signals using cross-correlation. The function returns the
delay in seconds. The function assumes that `data1` is the reference signal and `data2`
is the signal that is delayed with respect to `data1`. The funciton assumes that `data1`
and `data2` are sampled at the same rate and are of same length. `dt` is the time
resolution of the data.
"""
function find_delay(
    data1::Vector{Float64},
    data2::Vector{Float64},
    dt::Float64,
)::Float64
    corr = xcorr(data1, data2)
    return (argmax(abs.(corr)) - length(data1) + 1) * dt
end

"""
    gi_model(x::Vector{Float64}, p::Vector{Float64})::Vector{Float64}

Non linear gas injection model. The model is given by:

``\\Gamma = p_1 (\\sqrt{x^2 p_2^2 + 1} - 1)``

where `x` is the input voltage in Volts and ``\\Gamma`` is the flow rate in Pa m³/ s
"""
function gi_model(x::Vector{Float64}, p::Vector{Float64})::Vector{Float64}
    return p[1] .* (sqrt.((x .* p[2]) .^ 2 .+ 1) .- 1)  # Pa m³/ s
end

"""
    int_gi_model(x::Vector{Float64}, p::Vector{Float64})::Vector{Float64}

Cumulative sum of gas injection model. This function models the accumulated gas in the
vessel volume. Note that this is a numerical cumulative sum and does not change the
units of output. Returned output is in Pa m³/ s.
"""
function int_gi_model(x::Vector{Float64}, p::Vector{Float64})::Vector{Float64}
    ret_val = zeros(length(x))
    for ii ∈ eachindex(x)
        ret_val[ii] = sum(gi_model(x[1:ii], p))
    end
    return ret_val # Pa m³/ s
end

"""
    get_gas_injection_response(
        cmd::Vector{Float64},
        cmd_tt::Vector{Float64},
        P_ves::Vector{Float64},
        P_ves_tt::Vector{Float64},
        V_ves::Float64;
        resample_dt::Float64=0.02,
        gi_fit_guess::Vector{Float64}=[1.0, 1.0],
        response_curve_voltage::Vector{Float64}=collect(0:0.001:11),
    )::Tuple{IMAS.gas_injection__valve___response_curve, Dict{Symbol, Any}}

Function to fit a gas calibration shot data to a gas injection model. It requires 2 set
of data:

  - `cmd`, and `cmd_tt` are the command sent in Volts and the corresponding time axis.
  - `P_ves` and `P_ves_tt` are the pressure in the vessel and the corresponding time axis.

`V_ves` is the volume of the vessel. The function returns a tuple of `response_curve`
object and a valve model dictionary that can be used in
[`compute_gas_injection!`](@ref) or [`compute_gas_injection`](@ref) functions.

Optional keyword arguments:

  - `resample_dt`: Time resolution to resample the data.
  - `gi_fit_guess`: Initial guess for the gas injection model parameters.
  - `response_curve_voltage`: Voltage axis for storing the response curve.
"""
function get_gas_injection_response(
    cmd::Vector{Float64},
    cmd_tt::Vector{Float64},
    P_ves::Vector{Float64},
    P_ves_tt::Vector{Float64},
    V_ves::Float64;
    resample_dt::Float64=0.02,
    gi_fit_guess::Vector{Float64}=[1.0, 1.0],
    response_curve_voltage::Vector{Float64}=collect(0:0.001:11),
)::Tuple{IMAS.gas_injection__valve___response_curve, Dict{Symbol, Any}}
    tt =
        collect(min(cmd_tt[1], P_ves_tt[1]):resample_dt:max(cmd_tt[end], P_ves_tt[end]))
    cmd = downsample_smooth(cmd_tt, cmd, tt)
    P_ves = downsample_smooth(P_ves_tt, P_ves, tt)

    # Find latency between cmd and P_ves (derivative of P is easier to use for latency)
    dif_P_ves = [0; diff(P_ves)] / resample_dt # Pa / s
    latency = find_delay(dif_P_ves, cmd, resample_dt)
    latency_len = round(Int, latency / resample_dt)

    # Shift cmd forward by latency to match P_ves
    cmd = [zeros(latency_len); cmd[1:end-latency_len]]

    # Accumulated gas in vessel
    acc_gas = V_ves * P_ves  # m^3 * Pa

    # Fit the gas injection model to the data
    fit = curve_fit(
        int_gi_model,
        cmd[acc_gas.>0],
        acc_gas[acc_gas.>0] / resample_dt,
        gi_fit_guess,
    )

    # Calculate the response curve
    response_curve = IMAS.gas_injection__valve___response_curve()
    response_curve.voltage = response_curve_voltage
    response_curve.flow_rate = gi_model(response_curve_voltage, coef(fit))

    valve_model = Dict{Symbol, Any}(:latency => latency)

    return response_curve, valve_model
end

"""
    get_required_gas_cmd(
            required_flow_rate::Float64,
            response_curve::IMAS.gas_injection__valve___response_curve,
        )::Float64

Function to get the required gas command voltage to achieve the required flow rate.
"""
function get_required_gas_cmd(
    required_flow_rate::Float64,
    response_curve::IMAS.gas_injection__valve___response_curve,
)::Float64
    gas_cmd = linear_interpolation(response_curve.flow_rate, response_curve.voltage)
    return gas_cmd(required_flow_rate)
end

"""
    get_required_gas_cmd(
        required_flow_rate::Vector{Float64},
        response_curve::IMAS.gas_injection__valve___response_curve,
    )::Vector{Float64}

Function to get the required gas command voltage to achieve the required flow rate.
"""
function get_required_gas_cmd(
    required_flow_rate::Vector{Float64},
    response_curve::IMAS.gas_injection__valve___response_curve,
)::Vector{Float64}
    gas_cmd = linear_interpolation(response_curve.flow_rate, response_curve.voltage)
    return gas_cmd.(required_flow_rate)
end
