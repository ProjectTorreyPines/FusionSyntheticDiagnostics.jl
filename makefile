SHELL := /bin/zsh
help:
	@echo "Help Menu"
	@echo
	@echo "make env_with_cloned_repo (or make r): Creates a Julia environment with the cloned repositories"
	@echo "make env_with_git_url (or make u): Creates a Julia environment with the git urls without creating local clones"
	@echo "make clean: Deletes Project.toml and Manifest.toml for a fresh start"
	@echo

env_with_cloned_repo r:
	@echo "Pulling sample files using dvc"
	-dvc pull
	@echo "Creating Julia environment by creating local clones of dependent repositories"
	@echo "Cloning the repositories and generating Manifest.toml"
	-dn=$(shell dirname $(shell pwd)); \
	if [[ "$${dn:(-10)}" == ".julia/dev" ]]; then ext="" ; else ext=".jl";fi; \
	git clone "git@github.com:ProjectTorreyPines/OMAS.jl.git" ../OMAS$${ext}; \
	git clone "git@github.com:ProjectTorreyPines/GGDUtils.jl.git" ../GGDUtils$${ext}; \
	julia --project=. -e 'using Pkg; Pkg.rm(["OMAS", "GGDUtils"]); Pkg.develop(path="../OMAS'$${ext}'"); Pkg.develop(path="../GGDUtils'$${ext}'"); Pkg.instantiate()'

env_with_git_url u:
	@echo "Pulling sample files using dvc"
	-dvc pull
	@echo "Creating Julia environment with the git urls without creating local clones"
	@echo "Generating Project.toml and Manifest.toml"
	julia --project=. -e 'using Pkg; Pkg.rm(["OMAS", "GGDUtils"]); Pkg.add(url="git@github.com:ProjectTorreyPines/OMAS.jl.git", rev="master");  Pkg.add(url="git@github.com:ProjectTorreyPines/GGDUtils.jl.git", rev="master"); Pkg.instantiate()'

clean:
	@echo "Deleting Manifest.toml"
	-rm Manifest.toml

