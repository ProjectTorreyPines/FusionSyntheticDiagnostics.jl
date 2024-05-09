# SynthDiag.jl

![Format Check](https://github.com/ProjectTorreyPines/SynthDiag.jl/actions/workflows/format_check.yml/badge.svg)
![Docs](https://github.com/ProjectTorreyPines/SynthDiag.jl/actions/workflows/make_docs.yml/badge.svg)
![Tests](https://github.com/ProjectTorreyPines/SynthDiag.jl/actions/workflows/test.yml/badge.svg)

Package that defines synthetic diagnostics to generate sensor data based on plasma profile and syntehtic actuators. For installation and usage instructions, see the [online documentation](https://projecttorreypines.github.io/SynthDiag.jl/stable). For documentation on under development branch, see [dev online documentation](https://projecttorreypines.github.io/SynthDiag.jl/dev).

## Building julia environment (installation)

After cloning this repo, check the make menu:
```
SynthDiag.jl % make help
Help Menu

make env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories
make env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones
make clean: Deletes Project.toml and Manifest.toml for a fresh start
```

### make r
This option creates local copies of required private repositories at the same level as current repository and uses them in develop mode to create a Manifest.toml

### make u
This option uses url of required private repositories to create a static Manifest.toml attached to current master branches of these repositories.

### make clean
Deletes Manifest.toml so that environment can be recreated, to update or change the last used method.