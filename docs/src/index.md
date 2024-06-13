
# SynthDiag.jl 

```@contents
Pages = ["index.md"]
Depth = 5
```

## Installation

### Using make:
After cloning this repo, check the make menu:
```
SynthDiag.jl % make help
Help Menu

make env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories
make env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones
make clean: Deletes Project.toml and Manifest.toml for a fresh start
```

#### make r
This option creates local copies of required private repositories at the same level as current repository and uses them in develop mode to create a Manifest.toml

#### make u
This option uses url of required private repositories to create a static Manifest.toml attached to current master branches of these repositories.

#### make clean
Deletes Manifest.toml so that environment can be recreated, to update or change the last used method.

### Using Julia REPL and installing using Github url

Or, in julia REPL:
```julia
julia> using Pkg;
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/IMASDD.jl.git");
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/GGDUtils.jl.git");
julia> Pkg.add(; url="https://github.com/ProjectTorreyPines/SynthDiag.jl.git");
julia> Pkg.instantiate()
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
SynthDiag.FoV
SynthDiag.get_FoV
SynthDiag.get_line
SynthDiag.get_angle
SynthDiag.get_angle_bisector
SynthDiag.compute_intersection
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
