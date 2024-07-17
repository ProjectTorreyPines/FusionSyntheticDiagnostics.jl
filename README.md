# SynthDiag.jl

![Format Check](https://github.com/ProjectTorreyPines/SynthDiag.jl/actions/workflows/format_check.yml/badge.svg)
![Docs](https://github.com/ProjectTorreyPines/SynthDiag.jl/actions/workflows/make_docs.yml/badge.svg)
![Tests](https://github.com/ProjectTorreyPines/SynthDiag.jl/actions/workflows/test.yml/badge.svg)

Package that defines synthetic diagnostics to generate sensor data based on plasma profile and syntehtic actuators. For installation and usage instructions, see the [online documentation](https://projecttorreypines.github.io/SynthDiag.jl/stable). For documentation on under development branch, see [dev online documentation](https://projecttorreypines.github.io/SynthDiag.jl/dev).

## Installation

SynthDiag is registered with public repository [FuseRegistry](https://github.com/ProjectTorreyPines/FuseRegistry.jl/). For installation:

```
using Pkg
Pkg.Registry.add(RegistrySpec(url="https://github.com/ProjectTorreyPines/FuseRegistry.jl.git"))
Pkg.Registry.add("General")
Pkg.add("SynthDiag")
```