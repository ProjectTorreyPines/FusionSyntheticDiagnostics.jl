
# SynthDiag.jl 

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

SynthDiag is registered with public repository [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). First [install Julia](https://github.com/JuliaLang/juliaup?tab=readme-ov-file#juliaup---julia-version-manager), then:

```julia
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("SynthDiag")
```

## Synthetic Diagnostics

### Interferometer

Synthetic interferometer can be added using IMAS compatible `JSON` file that describes the metadata and the line of sight of chords. On computation, the interferometer uses edge profiles and core profiles data in the IMAS IDS to numerically integrate electron density along the line of sight and returns data in IMAS IDS interferometer object.

```@docs
add_interferometer!
compute_interferometer!
```

### Langmuir Probes

Langmuir probes can be added using IMAS compatible `JSON` file that describes the metadata and positions of embedded or reciprocating probes. Computation is currently supported for embedded langmuir probes only which uses the edge profiles data to report the edge electron and average ion temperature and electron density. If plasma potential and probe biasing information is available, it will use a langmuir probe current model to also calculate the ion saturaton current and current density as reported by a typical probe in IMAS data format.

```@docs
add_langmuir_probes!
compute_langmuir_probes!
langmuir_probe_current
```

## Synthetic Actuators

### Gas Injection

Gas valves can be added using IMAS compatible `JSON` file that describes the metadata, response curve, and positions of gas valves. If `ids.gas_injection.valve[:].voltage.data` is present, the gas output flow rate is calculated. Additionally, if a valve model dictionary is passed, more realistic actuation can be modeled that includes second order low pass effect, latency, and dribble effect.

```@docs
add_gas_injection!
compute_gas_injection
compute_gas_injection!
SynthDiag.get_lpf
SynthDiag.dribble
SynthDiag.downsample_smooth
SynthDiag.find_delay
get_gas_injection_response
SynthDiag.gi_model
SynthDiag.int_gi_model
get_required_gas_cmd
```

## Noise model

```@docs
Noise
generate_noise
generate_noise!
```
