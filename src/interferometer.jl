import PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e
import SD4SOLPS: fill_in_extrapolated_core_profile, check_rho_1d, add_rho_to_equilibrium

default_ifo = "$(@__DIR__)/default_interferometer.yml"

"""
    add_interferometer(
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
    config::String=default_ifo,

)::OMAS.dd

Add interferometer to IMAS structure using either a YAML or JSON file and compute the
line integrated electron density if not present
"""
function add_interferometer(
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
    config::String=default_ifo,
)::OMAS.dd
    if endswith(config, ".yml") || endswith(config, ".yaml")
        OMAS.yaml2imas(config, ids)
    elseif endswith(config, ".json")
        OMAS.json2imas(config, ids)
    end
    compute_interferometer(ids)
    return ids
end

"""
    add_interferometer(
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
    config::Dict{Symbol, Any},

)::OMAS.dd

Add interferometer to IMAS structure using a Dict and compute the line integrated
electron density if not present
"""
function add_interferometer(
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
    config::Dict{Symbol, Any},
)::OMAS.dd
    OMAS.dict2imas(config, ids)
    compute_interferometer(ids)
    return ids
end

"""
    compute_interferometer(@nospecialize(ids::OMAS.dd))

  - Calculate phase_to_n_e_line if not present for each wavelength
  - Compute the line integrated electron density if not present
"""
function compute_interferometer(@nospecialize(ids::OMAS.dd))
    for ch ∈ ids.interferometer.channel
        k = [2 * π / ch.wavelength[ii].value for ii ∈ 1:2]
        for i1 ∈ 1:2
            lam = ch.wavelength[i1]
            i2 = i1 % 2 + 1
            if lam.phase_to_n_e_line == 0
                # Taken from https://doi.org/10.1063/1.1138037
                lam.phase_to_n_e_line =
                    2 * m_e * ϵ_0 * c_0^2 / e^2 * (k[i1] * k[i2]^2) /
                    (k[i2]^2 - k[i1]^2)
            end
        end
        # Special case when measurement has not been made but edge profile data exists
        if length(ch.n_e_line.time) == 0 && length(ids.edge_profiles.ggd) > 0
            epggd = ids.edge_profiles.ggd
            cpp1d = ids.core_profiles.profiles_1d
            ch.n_e_line.time = ch.n_e_line.value = zeros(length(epggd))
            # Check if equilibrium data contains rho, if not add it
            if !check_rho_1d(ids; time_slice=1)
                add_rho_to_equilibrium!(ids)
            end
            # Check if core_profile is available, if not extrapolate from edge profile
            if length(cpp1d) == 0
                fill_in_extrapolated_core_profile(ids, "electrons.density")
            elseif length(cpp1d) != length(epggd)
                error(
                    "Number of edge profiles time slices does not match number of " \
                    "core profile time slices",
                )
            end
            fix_eq_time_idx = length(ids.equilibrium.time_slice) == 1
            for ii ∈ eachindex(epggd)
                ch.n_e_line.time[ii] = epggd[ii].time
                eq_time_idx = fix_eq_time_idx ? 1 : ii
                # WIP Integrate electron density along line of sight and add to value
                ch.n_e_line.value[ii]
            end
        end
    end
end
