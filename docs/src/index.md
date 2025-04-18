
# FusionSyntheticDiagnostics.jl 

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

First [install Julia](https://github.com/JuliaLang/juliaup?tab=readme-ov-file#juliaup---julia-version-manager), then:

```julia
using Pkg
Pkg.add("FusionSyntheticDiagnostics")
```

## Synthetic Diagnostics

### Bolometer

Synthetic bolometer can be added using IMAS compatible `JSON` file that describes the metadata and the line of sight of chords with information on apertures and detectors. On computation, the bolomter uses radiation data in the IMAS IDS to numerically integrate total radiation falling on each detector taking into account the apertures in the path.

```@docs
add_bolometer!
compute_bolometer!
```

Several useful geometrical functions are defined and used here.
```@docs
FusionSyntheticDiagnostics.FoV
FusionSyntheticDiagnostics.get_FoV
FusionSyntheticDiagnostics.get_line
FusionSyntheticDiagnostics.get_angle
FusionSyntheticDiagnostics.get_angle_bisector
FusionSyntheticDiagnostics.compute_intersection
```

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

### Magnetics

Magnetics can be added using IMAS compatible `JSON` file that describes the metadata and positions of Magntic field probes for poloidal and toroidal fields, flux loops, Rogowski coils, and shunts. Computation is currently supported for poloidal magnetic field probes and flux loops only which uses [IMAS physics flux-surfaces](https://projecttorreypines.github.io/IMAS.jl/dev/api/#Physics-flux-surfaces) functions to compute the fields and report them. As such, these diagnostic implementations are only reporting stored IMAS data in the format a real diagnostic would do, but they do not involve physical modeling or noise of the diagnostic as of now.

```@docs
add_magnetics!
compute_magnetics!
```

## Synthetic Actuators

### Gas Injection

Gas valves can be added using IMAS compatible `JSON` file that describes the metadata, response curve, and positions of gas valves. If `ids.gas_injection.valve[:].voltage.data` is present, the gas output flow rate is calculated. Additionally, if a valve model dictionary is passed, more realistic actuation can be modeled that includes second order low pass effect, latency, and dribble effect.

```@docs
add_gas_injection!
compute_gas_injection
compute_gas_injection!
FusionSyntheticDiagnostics.get_lpf
FusionSyntheticDiagnostics.dribble
FusionSyntheticDiagnostics.downsample_smooth
FusionSyntheticDiagnostics.find_delay
get_gas_injection_response
FusionSyntheticDiagnostics.gi_model
FusionSyntheticDiagnostics.int_gi_model
get_required_gas_cmd
```

## Noise model

```@docs
Noise
generate_noise
generate_noise!
```
