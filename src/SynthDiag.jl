module SynthDiag

using IMASDD: IMASDD
using StaticArrays
import PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e
import QuadGK: quadgk, BatchIntegrand
import GGDUtils:
    interp, get_grid_subset, get_subset_boundary, subset_do, get_TPS_mats
default_ifo = "$(@__DIR__)/default_interferometer.json"
default_gi = "$(@__DIR__)/default_gas_injection.json"

include("$(@__DIR__)/noise.jl")

include("$(@__DIR__)/interferometer.jl")

include("$(@__DIR__)/langmuir_probes.jl")

include("$(@__DIR__)/gas_injection.jl")

end # module SynthDiag
