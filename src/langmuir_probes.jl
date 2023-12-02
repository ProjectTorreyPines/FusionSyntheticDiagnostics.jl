import GGDUtils: interp, get_types_with

default_lp = "$(@__DIR__)/default_langmuir_probe.json"

"""
    add_langmuir_probes!(
    config::String=default_lp,
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),

)::OMAS.dd

Add langmuir_probes positions and other paramters to ids structure
"""
function add_langmuir_probes!(
    config::String=default_lp,
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
)::OMAS.dd
    if endswith(config, ".json")
        OMAS.json2imas(config, ids)
    else
        error("Only JSON files are supported.")
    end
    compute_langmuir_probes(ids)
    return ids
end

function compute_langmuir_probes(ids::OMAS.dd)
    epggd = ids.edge_profiles.ggd
    nt = length(epggd)
    fix_ep_grid_ggd_idx = length(ids.edge_profiles.grid_ggd) == 1
    ep_grid_ggd = ids.edge_profiles.grid_ggd[1]
    # Using -5 for now to use the SOLPS edge profile grid only
    TPS_mats_all_cells = get_TPS_mats(ep_grid_ggd, -5)
    # Get the edge profile interpolation functions
    ep_n_e = interp(epggd[1].electrons.density, TPS_mats_all_cells, 5)
    ep_t_e = interp(epggd[1].electrons.temperature, TPS_mats_all_cells, 5)

    # Initialize langmuir probe data
    init_data!.(ids.langmuir_probes.embedded, nt)

    for ii ∈ eachindex(epggd)
        if !fix_ep_grid_ggd_idx && ii > 1
            # If grid_ggd is evolving with time, update boundaries
            ep_grid_ggd = ids.edge_profiles.grid_ggd[ii]
            TPS_mats_all_cells = get_TPS_mats(ep_grid_ggd, -5)
            # Update the edge profile interpolation functions
            ep_n_e = interp(epggd[ii].electrons.density, TPS_mats_all_cells, 5)
            ep_t_e = interp(epggd[ii].electrons.temperature, TPS_mats_all_cells, 5)
        end

        for emb_lp ∈ ids.langmuir_probes.embedded
            emb_lp.time[ii] = epggd[ii].time
            emb_lp.n_e.data[ii] = ep_n_e(emb_lp.position.r, emb_lp.position.z)
            emb_lp.t_e.data[ii] = ep_t_e(emb_lp.position.r, emb_lp.position.z)
        end
    end
end

lp_data_types = Union{get_types_with(OMAS.langmuir_probes, :data)...}

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
        OMAS.langmuir_probes__embedded,
        OMAS.langmuir_probes__reciprocating___plunge,
        OMAS.langmuir_probes__reciprocating___plunge___collector,
    },
    nt::Int64,

)

Initialize each data field in an embedded langmuir probe.

For reciprocating probes, initialize each plunge needs to be initialized separately as
they can have different lengths.
"""
function init_data!(
    q::Union{
        OMAS.langmuir_probes__embedded,
        OMAS.langmuir_probes__reciprocating___plunge,
        OMAS.langmuir_probes__reciprocating___plunge___collector,
    },
    nt::Int64,
)
    for f ∈ fieldnames(typeof(q))
        fobj = getfield(q, f)
        type_fobj = typeof(fobj)
        if type_fobj <: lp_data_types
            init_data!(fobj, nt)
        elseif f == :time || f == :time_within_plunge
            setproperty!(q, f, zeros(Float64, nt))
        elseif type_fobj <: AbstractArray
            if eltype(fobj) <: OMAS.langmuir_probes__reciprocating___plunge___collector
                init_data!.(fobj, nt)
            end
        end
    end
end

"""
    init_data!(
    q::Union{
        OMAS.langmuir_probes__embedded,
        OMAS.langmuir_probes__reciprocating___plunge,
        OMAS.langmuir_probes__reciprocating___plunge___collector,
    },
    t::Vector{Float64},

)

Initialize probes data fields along with a time vector t for the time storing field
which is :time for embedded probes and :time_within_plunge for reciprocating probes.
"""
function init_data!(
    q::Union{
        OMAS.langmuir_probes__embedded,
        OMAS.langmuir_probes__reciprocating___plunge,
        OMAS.langmuir_probes__reciprocating___plunge___collector,
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
    init_data!(q::OMAS.langmuir_probes__reciprocating, nt::Int64)

Initialize each plunge in a reciprocating probe with nt length zeros.
"""
init_data!(q::OMAS.langmuir_probes__reciprocating, nt::Int64) =
    for p ∈ q.plunge
        init_data!(p, nt)
    end

"""
    init_data!(q::OMAS.langmuir_probes__reciprocating, t::Vector{Float64})

Initialize each plunge in a reciprocating probe with a time vector t.
"""
function init_data!(q::OMAS.langmuir_probes__reciprocating, t::Vector{Float64})
    for p ∈ q.plunge
        init_data!(p, t)
    end
end