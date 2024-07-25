using SynthDiag: add_interferometer!, add_langmuir_probes!, Noise,
    OverwriteAttemptError, add_gas_injection!, compute_gas_injection,
    get_gas_injection_response
using SynthDiag.IMAS: json2imas, dd
using DelimitedFiles: readdlm

println("-----------------------------------------------------------------------------")
ids = json2imas("$(@__DIR__)/../samples/time_dep_edge_profiles_with_equilibrium.json")
ids2 = json2imas("$(@__DIR__)/../samples/time_dep_edge_profiles_with_equilibrium.json")
print("add_interferometer!() time with compilation: ")
@time add_interferometer!(
    "$(@__DIR__)/../src/default_interferometer.json",
    ids;
    n_e_gsi=-5,
)
print("add_interferometer!() time (true runtime): ")
@time add_interferometer!(
    "$(@__DIR__)/../src/default_interferometer.json",
    ids2;
    n_e_gsi=-5,
)

println("-----------------------------------------------------------------------------")
# Assume a 5% noise level in ne values
ff = 0.0:0.1:1000
df = ff[2] - ff[1]
lpf = [f < 10.0 ? 1.0 : 100.0 * f^(-2) for f ∈ ff]
ne_noise_power = 1e14 * df * lpf
ne_noise = Noise(ne_noise_power, ff)
print("add_langmuir_probes!() time with compilation: ")
@time add_langmuir_probes!(
    "$(@__DIR__)/../src/default_langmuir_probes.json",
    ids;
    ne_noise=ne_noise,
    n_e_gsi=-5,
)
print("add_langmuir_probes!() time (true runtime): ")
@time add_langmuir_probes!(
    "$(@__DIR__)/../src/default_langmuir_probes.json",
    ids2;
    ne_noise=ne_noise,
    n_e_gsi=-5,
)

println("-----------------------------------------------------------------------------")
config = "$(@__DIR__)/../src/default_gas_injection.json"
excitation(t) = 0.6
ids = dd();
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

cmd_data = readdlm("$(@__DIR__)/../samples/gi193607_gacgasd.txt"; comments=true)
P_ves_data =
    readdlm("$(@__DIR__)/../samples/gi193607_pcm240tor.txt"; comments=true)
cmd_tt = cmd_data[:, 1] ./ 1000.0 # Convert to seconds
cmd = cmd_data[:, 2]
P_ves_tt = P_ves_data[:, 1] ./ 1000.0 # Convert to seconds
P_ves = P_ves_data[:, 2] * 1e-3 * 133.3223684211 # Convert to Pa
V_ves = 37.0 # m^3

print("get_gas_injection_response() time with compilation: ")
@time gasd_resp_curve, gasd_model =
    get_gas_injection_response(cmd, cmd_tt, P_ves, P_ves_tt, V_ves);
print("get_gas_injection_response() time (true runtime): ")
@time gasd_resp_curve, gasd_model =
    get_gas_injection_response(cmd, cmd_tt, P_ves, P_ves_tt, V_ves);

ids.gas_injection.valve[1].response_curve = gasd_resp_curve
# Adding made up time_constant and damping
gasd_model[:time_constant] = 0.3
gasd_model[:damping] = 0.8

valves = Dict("GASD" => gasd_model)

print("compute_gas_injection() time with compilation: ")
@time compute_gas_injection(ids; valves=valves)

ids2 = dd();
add_gas_injection!(config, ids2)

ids2.gas_injection.valve[1].voltage.time = tt
ids2.gas_injection.valve[1].voltage.data = zeros(nt)
ids2.gas_injection.valve[1].voltage.data[tstartind:tendind] .=
    excitation.(tt[tstartind:tendind])

ids2.gas_injection.valve[2].voltage.time = tt
ids2.gas_injection.valve[2].voltage.data = zeros(nt)
ids2.gas_injection.valve[2].voltage.data[tstartind2:tendind2] .=
    excitation.(tt[tstartind:tendind])

ids2.gas_injection.valve[1].response_curve = gasd_resp_curve

print("compute_gas_injection() time (true runtime): ")
@time compute_gas_injection(ids2; valves=valves)
