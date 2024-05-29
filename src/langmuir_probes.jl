import GGDUtils: interp, get_types_with
import PhysicalConstants.CODATA2018: m_e, m_u, e

export add_langmuir_probes!, compute_langmuir_probes!, langmuir_probe_current

default_lp = "$(@__DIR__)/default_langmuir_probes.json"

"""
    add_langmuir_probes!(
        config::Union{String, Dict{Symbol, Any}}=default_lp,
        @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
        overwrite=false, verbose=false, kwargs...,
    )::IMASDD.dd

Add langmuir probes positions and other parameters from `JSON` file or Julia `Dict` to
ids structure and compute langmuir probe outputs using edge profiles data. `kwargs` are
passed to [`compute_langmuir_probes!`](@ref).
"""
function add_langmuir_probes!(
    config::Union{String, Dict{Symbol, Any}}=default_lp,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite=false, verbose=false, kwargs...,
)::IMASDD.dd
    add_diagnostic!(
        config,
        :langmuir_probes,
        ids;
        overwrite=overwrite,
        verbose=verbose,
        channel=[:embedded, :reciprocating],
    )
    compute_langmuir_probes!(ids; kwargs...)
    return ids
end

"""
    compute_langmuir_probes!(
        ids::IMASDD.dd;
        v_plasma::Union{Float64, Vector{Float64}, Nothing}=nothing,
        v_floating::Union{Float64, Vector{Float64}, Nothing}=nothing,
        v_probe::Union{Float64, Vector{Float64}, Nothing}=nothing,
        v_plasma_noise::Union{Noise, Nothing}=nothing,
        Z_eff::Float64=1.0,
        m_i_amu::Float64=2.014, # Deuterium
        probe_type::Symbol=:cylindrical,
        ne_noise::Union{Noise, Nothing}=nothing,
        te_noise::Union{Noise, Nothing}=nothing,
        ti_noise::Union{Noise, Nothing}=nothing,
        n_e_gsi::Int=5,
    )

Compute langmuir probe outputs for all the embedded langmuir probes in the ids structure
using edge profiles data. If `v_plasma` or `v_floating` are provided along with
`v_probe`, then [`langmuir_probe_current()`](@ref) is used to compute the ion saturation
current and parallel current density and stored in the IDS. Noise models can be added
for electron density and temperature and ion average temperature as power spectral
densities using data type [`Noise`](@ref).

Input arguments:

  - `ids`: IMASDD.dd object
  - `v_plasma`: Plasma potential (in V) or a vector of plasma potentials for each time
    step
  - `v_floating`: Floating potential (in V) or a vector of floating potentials for each
    time step
  - `v_probe`: Probe potential (in V) or a vector of probe potentials for each time step
  - `v_plasma_noise`: Noise object for plasma potential
  - `Z_eff`: Effective charge of the ion
  - `m_i_amu`: Ion mass in atomic mass units (amu)
  - `probe_type`: Type of the probe, either :cylindrical or :spherical
  - `ne_noise`: Noise object for electron density
  - `te_noise`: Noise object for electron temperature
  - `ti_noise`: Noise object for average ion temperature
  - `n_e_gsi`: Grid subset index that is used in edge profiles data for electron density,
    electron temperature, and average ion temperature.
"""
function compute_langmuir_probes!(
    ids::IMASDD.dd;
    v_plasma::Union{Float64, Vector{Float64}, Nothing}=nothing,
    v_floating::Union{Float64, Vector{Float64}, Nothing}=nothing,
    v_probe::Union{Float64, Vector{Float64}, Nothing}=nothing,
    v_plasma_noise::Union{Noise, Nothing}=nothing,
    Z_eff::Float64=1.0,
    m_i_amu::Float64=2.014, # Deuterium
    probe_type::Symbol=:cylindrical,
    ne_noise::Union{Noise, Nothing}=nothing,
    te_noise::Union{Noise, Nothing}=nothing,
    ti_noise::Union{Noise, Nothing}=nothing,
    n_e_gsi::Int=5,
)
    epggd = ids.edge_profiles.ggd
    nt = length(epggd)
    fix_ep_grid_ggd_idx = length(ids.edge_profiles.grid_ggd) == 1
    ep_grid_ggd = ids.edge_profiles.grid_ggd[1]

    TPS_mats = get_TPS_mats(ep_grid_ggd, n_e_gsi)

    ep_n_e_list = Array{Function}(undef, nt)
    ep_t_e_list = Array{Function}(undef, nt)
    ep_t_i_list = Array{Function}(undef, nt)
    for ii ∈ eachindex(epggd)
        TPS_mats_ii = update_TPS_mats(ii, fix_ep_grid_ggd_idx, ids, n_e_gsi, TPS_mats)
        ep_n_e_list[ii] = interp(epggd[ii].electrons.density, TPS_mats_ii, n_e_gsi)
        ep_t_e_list[ii] = interp(epggd[ii].electrons.temperature, TPS_mats_ii, n_e_gsi)
        ep_t_i_list[ii] = interp(epggd[ii].t_i_average, TPS_mats_ii, n_e_gsi)
    end

    # Initialize langmuir probe data
    init_data!.(ids.langmuir_probes.embedded, nt)

    for ii ∈ eachindex(epggd)
        ep_n_e = ep_n_e_list[ii]
        ep_t_e = ep_t_e_list[ii]
        ep_t_i = ep_t_i_list[ii]

        for emb_lp ∈ ids.langmuir_probes.embedded
            emb_lp.time[ii] = epggd[ii].time
            emb_lp.n_e.data[ii] = ep_n_e(emb_lp.position.r, emb_lp.position.z)
            emb_lp.t_e.data[ii] = ep_t_e(emb_lp.position.r, emb_lp.position.z)
            emb_lp.t_i.data[ii] = ep_t_i(emb_lp.position.r, emb_lp.position.z)
        end
    end
    if !isnothing(ne_noise)
        for emb_lp ∈ ids.langmuir_probes.embedded
            emb_lp.n_e.data .+= generate_noise(ne_noise, emb_lp.time)
        end
    end
    if !isnothing(te_noise)
        for emb_lp ∈ ids.langmuir_probes.embedded
            emb_lp.t_e.data .+= generate_noise(te_noise, emb_lp.time)
        end
    end
    if !isnothing(ti_noise)
        for emb_lp ∈ ids.langmuir_probes.embedded
            emb_lp.t_i.data .+= generate_noise(ti_noise, emb_lp.time)
        end
    end

    if !isnothing(v_plasma)
        for emb_lp ∈ ids.langmuir_probes.embedded
            if length(v_plasma) == 1
                emb_lp.v_plasma.data .= v_plasma * ones(nt)
            elseif length(v_plasma) == nt
                emb_lp.v_plasma.data .= v_plasma
            else
                error(
                    "v_plasma must be a scalar or a vector of length of size of " *
                    " edge_profiles.ggd corresponding to edge_profiles.ggd[:].time",
                )
            end
        end
        if !isnothing(v_plasma_noise)
            for emb_lp ∈ ids.langmuir_probes.embedded
                emb_lp.v_plasma.data .+= generate_noise(v_plasma_noise, emb_lp.time)
            end
        end
    end
    if !isnothing(v_floating)
        for emb_lp ∈ ids.langmuir_probes.embedded
            if length(v_floating) == 1
                emb_lp.v_floating.data .= v_floating * ones(nt)
            elseif length(v_floating) == nt
                emb_lp.v_floating.data .= v_floating
            else
                error(
                    "v_floating must be a scalar or a vector of length of size of " *
                    " edge_profiles.ggd corresponding to edge_profiles.ggd[:].time",
                )
            end
        end
    end

    for emb_lp ∈ ids.langmuir_probes.embedded
        if !isnothing(v_plasma) && isnothing(v_floating)
            emb_lp.v_floating.data = Te .* log.(abs.(I_isc ./ I_esc))
        end
        if !isnothing(v_probe) && !IMASDD.ismissing(emb_lp.v_plasma.data)
            i_probe =
                langmuir_probe_current.(
                    v_probe - emb_lp.v_plasma.data,
                    emb_lp.t_e.data,
                    emb_lp.t_i.data,
                    emb_lp.n_e.data,
                    emb_lp.surface_area,
                    m_i_amu,
                    Z_eff;
                    probe_type=probe_type,
                )
        else
            i_probe =
                langmuir_probe_current.(
                    -200.0,                        # Assume 200 V negative bias
                    emb_lp.t_e.data,
                    emb_lp.t_i.data,
                    emb_lp.n_e.data,
                    emb_lp.surface_area,
                    m_i_amu,
                    Z_eff;
                    probe_type=probe_type,
                )
        end
        emb_lp.ion_saturation_current.data = i_probe
        emb_lp.j_i_parallel.data = i_probe ./ emb_lp.surface_area
    end
end

"""
    langmuir_probe_current(
        Δv::Float64,
        Te::Float64,
        Ti=::Float64,
        ne::Float64,
        A::Float64,
        m_i_amu=2.04,
        Z_eff::Float64=1.0;
        probe_type::Symbol=:cylindrical,
    )::Float64

Langmuir Probe current model for cylindrical and spherical probes.

For `Δv` <= 0:

``I_{probe} = -I_{isc} (1 - \\frac{Δv}{T_i})^{ex} - I_{esc} e^{\\frac{Δv}{T_e}}``

For `Δv` > 0:

``I_{probe} = -I_{esc} (1 + \\frac{Δv}{T_e})^{ex} - I_{isc} e^{-\\frac{Δv}{T_i}}``

where

``I_{isc} = \\frac{1}{2} A q_i n_e \\sqrt{e (T_e + T_i) / m_i}``

``I_{esc} = -\\frac{1}{4} A e n_e \\sqrt{8 e T_e / m_e}``

where `A` is the effective surface area of the probe, `q_i`(`Z_eff`e) is the total
charge of the ion, and `ex` is 1 for spherical and 0.5 for cylindrical probes.

Input arguments:

  - `Δv`: Bias voltage (`v_probe` - `v_plasma`)
  - `Te`: Electron temperature (in eV)
  - `Ti`: Ion temperature (in eV)
  - `ne`: Electron density (in ``m^{-3}``)
  - `A`: Effective surface area of the probe (in ``m^{2}``)
  - `m_i_amu`: Ion mass in atomic mass units (amu)
  - `Z_eff`: Effective charge of the ion (in units of elementary charge e)
  - `probe_type`: Type of the probe, either :cylindrical or :spherical

Ref: L. Conde, "[An introduction to Langmuir probe diagnostics of plasmas](https://api.semanticscholar.org/CorpusID:53622081)" (2011)
"""
function langmuir_probe_current(
    Δv::Float64,
    Te::Float64,
    Ti::Float64,
    ne::Float64,
    A::Float64,
    m_i_amu::Float64=2.04,
    Z_eff::Float64=1.0;
    probe_type::Symbol=:cylindrical,
)::Float64
    m_i = m_i_amu * m_u.val
    q_i = Z_eff * e.val
    I_isc = 0.5 * A * q_i * ne .* sqrt.(e.val .* (Te .+ Ti) ./ m_i)
    I_esc = -0.25 * A * e.val * ne .* sqrt.(8 * e.val .* Te ./ m_e.val)
    if probe_type == :spherical
        ex = 1
    elseif probe_type == :cylindrical
        ex = 0.5
    else
        error("Probe_type must be either :cylindrical or :spherical")
    end
    if Δv <= 0
        i_probe = -I_isc * (1 - Δv / Ti)^ex - I_esc * exp(Δv / Te)
    else
        i_probe = -I_esc * (1 + Δv / Te)^ex - I_isc * exp(-Δv / Ti)
    end
    return i_probe
end

lp_data_types = Union{get_types_with(IMASDD.langmuir_probes, :data)...}

"""
    init_data!(q::lp_data_types, nt::Int64)

Initialize langmuir probe data field for a measurement length of nt
"""
function init_data!(q::lp_data_types, nt::Int64)
    q.data = zeros(Float64, nt)
    q.validity_timed = zeros(Int64, nt)
    return q.validity = 0
end

"""
    init_data!(
        q::Union{
            IMASDD.langmuir_probes__embedded,
            IMASDD.langmuir_probes__reciprocating___plunge,
            IMASDD.langmuir_probes__reciprocating___plunge___collector,
        },
        nt::Int64,
    )

Initialize each data field in an embedded langmuir probe.

For reciprocating probes, initialize each plunge needs to be initialized separately as
they can have different lengths.
"""
function init_data!(
    q::Union{
        IMASDD.langmuir_probes__embedded,
        IMASDD.langmuir_probes__reciprocating___plunge,
        IMASDD.langmuir_probes__reciprocating___plunge___collector,
    },
    nt::Int64,
)
    # make sure to initialize time before any other field
    for f ∈ [:time; collect(fieldnames(typeof(q)))]
        fobj = getfield(q, f)
        type_fobj = typeof(fobj)
        if type_fobj <: lp_data_types
            init_data!(fobj, nt)
        elseif f == :time || f == :time_within_plunge
            setproperty!(q, f, zeros(Float64, nt))
        elseif type_fobj <: AbstractArray
            if eltype(fobj) <:
               IMASDD.langmuir_probes__reciprocating___plunge___collector
                init_data!.(fobj, nt)
            end
        end
    end
end

"""
    init_data!(
        q::Union{
            IMASDD.langmuir_probes__embedded,
            IMASDD.langmuir_probes__reciprocating___plunge,
            IMASDD.langmuir_probes__reciprocating___plunge___collector,
        },
        t::Vector{Float64},
    )

Initialize probes data fields along with a time vector t for the time storing field
which is :time for embedded probes and :time_within_plunge for reciprocating probes.
"""
function init_data!(
    q::Union{
        IMASDD.langmuir_probes__embedded,
        IMASDD.langmuir_probes__reciprocating___plunge,
        IMASDD.langmuir_probes__reciprocating___plunge___collector,
    },
    t::Vector{Float64},
)
    nt = length(t)
    init_data!(q, nt)
    if :time ∈ fieldnames(typeof(q))
        q.time = t
    elseif :time_within_plunge ∈ fieldnames(typeof(q))
        q.time_within_plunge = t
    end
end

"""
    init_data!(q::IMASDD.langmuir_probes__reciprocating, nt::Int64)

Initialize each plunge in a reciprocating probe with nt length zeros.
"""
init_data!(q::IMASDD.langmuir_probes__reciprocating, nt::Int64) =
    for p ∈ q.plunge
        init_data!(p, nt)
    end

"""
    init_data!(q::IMASDD.langmuir_probes__reciprocating, t::Vector{Float64})

Initialize each plunge in a reciprocating probe with a time vector t.
"""
function init_data!(q::IMASDD.langmuir_probes__reciprocating, t::Vector{Float64})
    for p ∈ q.plunge
        init_data!(p, t)
    end
end
