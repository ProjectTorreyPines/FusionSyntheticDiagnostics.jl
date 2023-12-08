module SynthDiag

using OMAS: OMAS
using StaticArrays
import PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e
import QuadGK: quadgk, BatchIntegrand
import GGDUtils:
    interp, get_grid_subset_with_index, get_subset_boundary, subset_do, get_TPS_mats, is_inside

default_ifo = "$(@__DIR__)/default_interferometer.json"


include("$(@__DIR__)/noise.jl")

include("$(@__DIR__)/interferometer.jl")

include("$(@__DIR__)/langmuir_probes.jl")

end # module SynthDiag
