using SynthDiag: add_interferometer!, add_langmuir_probes!, Noise, OverwriteAttemptError
using SynthDiag.IMASDD: json2imas

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
