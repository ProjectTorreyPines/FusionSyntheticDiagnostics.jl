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
