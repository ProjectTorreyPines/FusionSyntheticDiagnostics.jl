struct OverwriteAttemptError <: Exception
    var::String
end

Base.showerror(io::IO, e::OverwriteAttemptError) = print(io, e.var)

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
