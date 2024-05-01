import GGDUtils: interp, get_types_with
using PhysicalConstants.CODATA2018: m_e, m_u, e

"""
    add_langmuir_probes!(
    config::String=default_lp,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite=false, verbose=false, kwargs...,

)::IMASDD.dd

Add langmuir_probes positions and other parameters from json file to ids structure
"""
function add_langmuir_probes!(
    config::String=default_lp,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite=false, verbose=false, kwargs...,
)::IMASDD.dd
    if endswith(config, ".json")
        config_dict = convert_strings_to_symbols(IMASDD.JSON.parsefile(config))
        add_langmuir_probes!(
            config_dict,
            ids;
            overwrite=overwrite,
            verbose=verbose,
            kwargs...,
        )
    else
        error("Only JSON files are supported.")
    end
    return ids
end

"""
    add_langmuir_probes!(
    config::Dict{Symbol, Any},
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite=false, verbose=false, kwargs...,

)::IMASDD.dd

Add langmuir_probes positions and other parameters from Dictionary to ids structure
"""
function add_langmuir_probes!(
    config::Dict{Symbol, Any},
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite=false, verbose=false, kwargs...,
)::IMASDD.dd
    # Check for duplicates
    if length(ids.langmuir_probes.embedded) > 0
        new_lps = Dict()
        if haskey(config, :langmuir_probes)
            if haskey(config[:langmuir_probes], :embedded)
                for emp_lps ∈ config[:langmuir_probes][:embedded]
                    new_lps[emp_lps[:name]] = emp_lps[:identifier]
                end
            end
            if haskey(config[:langmuir_probes], :reciprocating)
                for rec_lp ∈ config[:langmuir_probes][:reciprocating]
                    new_lps[rec_lp[:name]] = rec_lp[:identifier]
                end
            end
        else
            warning("Config does not have langmuir_probes in it. Skipping.")
            return ids
        end
        dup_inds = Dict(:emb_lp => [], :rec_lp => [])
        for (ii, emb_lp) ∈ enumerate(ids.langmuir_probes.embedded)
            if emb_lp.name in keys(new_lps) ||
               emb_lp.identifier in values(new_lps)
                append!(dup_inds[:emb_lp], ii)
            end
        end
        for (ii, rec_lp) ∈ enumerate(ids.langmuir_probes.reciprocating)
            if rec_lp.name in keys(new_lps) ||
               rec_lp.identifier in values(new_lps)
                append!(dup_inds[:rec_lp], ii)
            end
        end
        if overwrite
            for ii ∈ reverse(dup_inds[:emb_lp])
                println(
                    "Overwriting embedded langmuir_probe ",
                    "$(ids.langmuir_probes.embedded[ii].name)...",
                )
                deleteat!(ids.langmuir_probes.embedded, ii)
            end
            for ii ∈ reverse(dup_inds[:rec_lp])
                println(
                    "Overwriting reciprocating langmuir_probes ",
                    "$(ids.langmuir_probes.reciprocating[ii].name)...",
                )
                deleteat!(ids.langmuir_probes.reciprocating, ii)
            end
        else
            if length(dup_inds[:emb_lp]) + length(dup_inds[:rec_lp]) > 0
                err_msg =
                    "Duplicate langmuir_probes embeddeds found with " *
                    "overlapping names or identifiers.\n" * "Identifier: Name\n"
                for ii ∈ dup_inds[:emb_lp]
                    err_msg *=
                        "$(ids.langmuir_probes.embedded[ii].identifier): " *
                        "$(ids.langmuir_probes.embedded[ii].name)\n"
                end
                for ii ∈ dup_inds[:rec_lp]
                    err_msg *=
                        "$(ids.langmuir_probes.reciprocating[ii].identifier): " *
                        "$(ids.langmuir_probes.reciprocating[ii].name)\n"
                end
                err_msg *= "Use overwrite=true to replace them."
                throw(OverwriteAttemptError(err_msg))
            end
        end
        config[:langmuir_probes] =
            mergewith(
                append!,
                IMASDD.imas2dict(ids.langmuir_probes),
                config[:langmuir_probes],
            )
    end
    IMASDD.dict2imas(config, ids; verbose=verbose)
    compute_langmuir_probes!(ids; kwargs...)
    return ids
end

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
    n_e_gsi::Int=5,
)
    epggd = ids.edge_profiles.ggd
    nt = length(epggd)
    fix_ep_grid_ggd_idx = length(ids.edge_profiles.grid_ggd) == 1
    ep_grid_ggd = ids.edge_profiles.grid_ggd[1]
    # TPS_mats_all_cells = get_TPS_mats(ep_grid_ggd, n_e_gsi)

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

    Input arguments:
    - Δv: Bias voltage (v_probe - v_plasma)
    - Te: Electron temperature (in eV)
    - Ti: Ion temperature (in eV)
    - ne: Electron density (in m^-3)
    - A: Surface area of the probe (in m^2)
    - m_i_amu: Ion mass in atomic mass units (amu)
    - Z_eff: Effective charge of the ion (in units of elementary charge e)
    - probe_type: Type of the probe, either :cylindrical or :spherical
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
