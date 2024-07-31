using SD4SOLPS: SD4SOLPS
using IMAS: IMAS
using Plots: Plots

eqdsk_file = "geqdsk_iter_small_sample"
sample_paths = [
    splitdir(pathof(SD4SOLPS))[1] * "/../sample",
]
core_method = "simple"
edge_method = "simple"
filename = splitdir(pathof(SD4SOLPS))[1] * "/../sd_input_data"
output_format = "json"
dd = SD4SOLPS.preparation(
    eqdsk_file,
    sample_paths;
    core_method=core_method,
    filename=filename,
    output_format=output_format,
    allow_boundary_flux_correction=true,
)

Rwall, Zwall, rwall, zwall, s, SOL, r, q = IMAS.mesher_HF(dd)
IMAS.particle_HF(dd.equilibrium.time_slice[1], SOL, rwall, zwall, r, q)
Plots.plot(s, Qwall)
