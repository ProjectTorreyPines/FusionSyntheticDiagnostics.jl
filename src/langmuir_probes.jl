import GGDUtils: interp, get_types_with
using PhysicalConstants.CODATA2018: m_e, m_u, e

"""
    add_langmuir_probes!(
    config::String=default_lp,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd(),

)::IMASDD.dd

Add langmuir_probes positions and other paramters to ids structure
"""
function add_langmuir_probes!(
    config::String=default_lp,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd(); kwargs...,
)::IMASDD.dd
    if endswith(config, ".json")
        IMASDD.json2imas(config, ids)
    else
        error("Only JSON files are supported.")
    end
    compute_langmuir_probes(ids; kwargs...)
    return ids
end

function compute_langmuir_probes(
    ids::IMASDD.dd;
    v_plasma::Union{Float64, Vector{Float64}, Nothing}=nothing,
    v_floating::Union{Float64, Vector{Float64}, Nothing}=nothing,
    v_probe::Union{Float64, Vector{Float64}, Nothing}=nothing,
    v_plasma_noise::Union{Noise, Nothing}=nothing,
    Z_eff::Float64=1.0,
    m_i_amu::Float64=2.014, # Deuterium
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

    # Get the edge profile interpolation functions
    # ep_n_e = interp(epggd[1].electrons.density, TPS_mats_all_cells, n_e_gsi)
    # ep_t_e = interp(epggd[1].electrons.temperature, TPS_mats_all_cells, n_e_gsi)
    # ep_t_i = interp(epggd[1].t_i_average, TPS_mats_all_cells, n_e_gsi)

    # Initialize langmuir probe data
    init_data!.(ids.langmuir_probes.embedded, nt)

    for ii ∈ eachindex(epggd)
        # if !fix_ep_grid_ggd_idx && ii > 1
        #     # If grid_ggd is evolving with time, update boundaries
        #     ep_grid_ggd = ids.edge_profiles.grid_ggd[ii]
        #     TPS_mats_all_cells = get_TPS_mats(ep_grid_ggd, n_e_gsi)
        #     # Update the edge profile interpolation functions
        #     ep_n_e = interp(epggd[ii].electrons.density, TPS_mats_all_cells, n_e_gsi)
        #     ep_t_e =
        #         interp(epggd[ii].electrons.temperature, TPS_mats_all_cells, n_e_gsi)
        # end

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

    m_i = m_i_amu * m_u.val
    q_i = Z_eff * e.val
    for emb_lp ∈ ids.langmuir_probes.embedded
        A = emb_lp.surface_area
        ne = emb_lp.n_e.data
        Te = emb_lp.t_e.data
        Ti = emb_lp.t_i.data
        I_isc = 0.5 * A * q_i * ne .* sqrt.(e.val .* (Te .+ Ti) ./ m_i)
        I_esc = -0.25 * A * e.val * ne .* sqrt.(8 * e.val .* Te ./ m_e)
        if !isnothing(v_plasma) && isnothing(v_floating)
            emb_lp.v_floating.data = Te .* log.(abs.(I_isc ./ I_esc))
        end
        if isnothing(v_probe) && !IMASDD.ismissing(emb_lp.v_plasma.data)
            emb_lp.ion_saturation_current.data = I_isc
            emb_lp.j_i_parallel.data = I_isc ./ A
        else
            if length(v_probe) == 1
                v_probe .= v_probe * ones(nt)
            elseif length(v_probe) != nt
                error(
                    "v_probe must be a scalar or a vector of length of size of " *
                    " edge_profiles.ggd corresponding to edge_profiles.ggd[:].time",
                )
            end
            # This is for cylindrical probe
            v_plasma = emb_lp.v_plasma.data
            i_probe =
                -I_isc .* sqrt.(1 .- (v_probe - v_plasma) ./ Ti) .+
                I_esc * exp.((v_probe .- v_plasma) ./ Te)
            emb_lp.ion_saturation_current.data = i_probe
            emb_lp.j_i_parallel.data = i_probe ./ A
        end
    end
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
