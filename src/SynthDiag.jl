module SynthDiag

using OMAS: OMAS

include("$(@__DIR__)/interferometer.jl")

include("$(@__DIR__)/langmuir_probes.jl")

end # module SynthDiag
