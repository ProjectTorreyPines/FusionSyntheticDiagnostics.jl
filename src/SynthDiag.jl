module SynthDiag

using IMASDD: IMASDD

include("$(@__DIR__)/noise.jl")

include("$(@__DIR__)/utils.jl")

include("$(@__DIR__)/interferometer.jl")

include("$(@__DIR__)/langmuir_probes.jl")

include("$(@__DIR__)/gas_injection.jl")

include("$(@__DIR__)/derived.jl")

include("$(@__DIR__)/magic.jl")

end # module SynthDiag
