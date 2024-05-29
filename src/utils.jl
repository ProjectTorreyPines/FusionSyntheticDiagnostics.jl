struct OverwriteAttemptError <: Exception
    var::String
end

Base.showerror(io::IO, e::OverwriteAttemptError) = print(io, e.var)

function add_diagnostic!(
    config::String,
    diagnostic::Symbol,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    channel::Union{Symbol, Vector{Symbol}}=:channel,
    overwrite=false,
    verbose=false,
)
    if endswith(config, ".json")
        config_dict = convert_strings_to_symbols(IMASDD.JSON.parsefile(config))
        add_diagnostic!(
            config_dict,
            diagnostic,
            ids;
            channel=channel,
            overwrite=overwrite,
            verbose=verbose,
        )
    else
        error("Only JSON files are supported.")
    end
    return ids
end

function add_diagnostic!(
    config::Dict{Symbol, Any},
    diagnostic::Symbol,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    channel::Union{Symbol, Vector{Symbol}}=:channel,
    overwrite=false,
    verbose=false,
)
    if isa(channel, Vector)
        for chan ∈ channel
            add_diagnostic!(
                config,
                diagnostic,
                ids;
                channel=chan,
                overwrite=overwrite,
                verbose=verbose,
            )
        end
        return ids
    end
    diag_name = String(diagnostic)
    if !haskey(config, diagnostic)
        warning("Config does not have " * diag_name * " in it. Skipping.")
        return ids
    end
    ids_diag = getproperty(ids, diagnostic)
    ids_diag_channels = getproperty(ids_diag, channel)
    if length(ids_diag_channels) > 0
        duplicate_indices = []
        new_channels = Dict(
            ch[:name] => ch[:identifier] for
            ch ∈ config[diagnostic][channel]
        )
        for (ii, ch) ∈ enumerate(ids_diag_channels)
            if ch.name in keys(new_channels) ||
               ch.identifier in values(new_channels)
                append!(duplicate_indices, ii)
            end
        end
        if overwrite
            for ii ∈ reverse(duplicate_indices)
                println(
                    "Overwriting " * diag_name * " channel ",
                    "$(ids_diag_channels[ii].name)...",
                )
                deleteat!(ids_diag_channels, ii)
            end
        else
            if length(duplicate_indices) > 0
                err_msg =
                    "Duplicate " * diag_name * " channels found with " *
                    "overlapping names or identifiers.\n" * "Identifier: Name\n"
                for ii ∈ duplicate_indices
                    err_msg *=
                        "$(ids_diag_channels[ii].identifier): " *
                        "$(ids_diag_channels[ii].name)\n"
                end
                err_msg *= "Use overwrite=true to replace them."
                throw(OverwriteAttemptError(err_msg))
            end
        end
        config[diagnostic] = mergewith(
            append!,
            IMASDD.imas2dict(ids_diag),
            config[diagnostic],
        )
    end
    IMASDD.dict2imas(config, ids; verbose=verbose)
    return ids
end

function convert_strings_to_symbols(d::Dict{String, Any})
    new_d = Dict{Symbol, Any}()
    for (k, v) ∈ d
        if isa(v, Dict{String, Any})
            new_d[Symbol(k)] = convert_strings_to_symbols(v)
        elseif isa(v, Vector{Any})
            if length(v) > 0
                if isa(v[1], Dict{String, Any})
                    new_d[Symbol(k)] = Vector{Dict{Symbol, Any}}(undef, length(v))
                    for ii ∈ eachindex(v)
                        new_d[Symbol(k)][ii] = convert_strings_to_symbols(v[ii])
                    end
                else
                    new_d[Symbol(k)] = v
                end
            else
                new_d[Symbol(k)] = v
            end
        else
            new_d[Symbol(k)] = v
        end
    end
    return new_d
end

@inline function xyz2rz(x::Float64, y::Float64, z::Float64)
    r = sqrt(x^2 + y^2)
    return r, z
end

function update_TPS_mats(ii, fix_ep_grid_ggd_idx, ids, gsi, TPS_mats)
    if !fix_ep_grid_ggd_idx
        ep_grid_ggd = ids.edge_profiles.grid_ggd[ii]
        return get_TPS_mats(ep_grid_ggd, gsi)
    else
        return TPS_mats
    end
end
