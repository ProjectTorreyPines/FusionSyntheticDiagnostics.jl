module FusionSyntheticDiagnostics

import IMASdd as IMAS

include("$(@__DIR__)/noise.jl")

include("$(@__DIR__)/utils.jl")

include("$(@__DIR__)/interferometer.jl")

include("$(@__DIR__)/langmuir_probes.jl")

include("$(@__DIR__)/gas_injection.jl")

end # module FusionSyntheticDiagnostics
