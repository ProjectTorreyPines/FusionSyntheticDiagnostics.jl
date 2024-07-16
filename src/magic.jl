# File for holding "magic" diagnostics: measurements and calculations based on
# information that probably wouldn't be accessible in a real plasma. Intended to
# give a reference to the "truth" (as much as we can call it that with our
# imperfect physics models) that we try to approach with the more realistic
# synthetic diagnostics.

# Some of these quantities could have counterparts estimated from realistic
# synthetic diagnostic data and must be distinguished. So, prefixing function
# names with magic_ is encouraged.

using GGDUtils:
    get_prop_with_grid_subset_index, get_subset_centers, get_grid_subset, interp

export magic_nesep

"""
    magic_nesep(
        dd::IMASDD.dd;
        grid_ggd_idx::Int64=1,
        space_number::Int64=1,
        cell_grid_subset::Int64=5,
        midplane_grid_subset::Int64=11,
        separatrix_grid_subset::Int64=16,
    )::Array{Float64}

Returns the value of electron density at the outboard midplane as a function of time.
The value is interpolated from the model.

Inputs:

  - dd: a data dictionary instance from IMASDD
  - grid_ggd_idx: index of the grid_ggd to use. Many cases will have only one grid_ggd.
  - space_number: space index to use. You may only have one space.
  - cell_grid_subset: index of the grid subset for cells. This is normally 5, but could be
    changed to, for example, -5 to index only cells that are from the original B2 mesh in
    the case that the mesh has been extended.
  - midplane_grid_subset: index of the grid subset for y-aligned faces at the outboard
    midplane. This is normally 11 and should not deviate from the standard without a reason.
  - separatrix_grid_subset: index of the subset for x-aligned faces at the separatrix. Leave
    this alone at 16 unless there's a reason for changing it.

Output:
nesep: an array of n_e,sep values in m^-3 vs time. If the dd has only one timeslice,
the array will only have one element but it will still be an array.
"""
function magic_nesep(
    dd::IMASDD.dd;
    grid_ggd_idx::Int64=1,
    space_number::Int64=1,
    cell_grid_subset::Int64=5,
    midplane_grid_subset::Int64=11,
    separatrix_grid_subset::Int64=16,
)::Array{Float64}
    # Form empty output array
    nslices = length(dd.edge_profiles.ggd)
    nesep = Array{Float64}(undef, nslices)

    # Get out references to basic ggd objects
    fix_ep_grid_ggd_idx = length(dd.edge_profiles.grid_ggd) == 1
    grid_ggd = dd.edge_profiles.grid_ggd[grid_ggd_idx]
    space = grid_ggd.space[space_number]

    rintersect, zintersect = find_OMP_RZ(
        grid_ggd;
        space_number=space_number,
        midplane_grid_subset=midplane_grid_subset,
        separatrix_grid_subset=separatrix_grid_subset,
    )

    tps_mats = get_TPS_mats(grid_ggd, cell_grid_subset)

    for it ∈ 1:nslices
        ggd = dd.edge_profiles.ggd[it]

        nesep[it] = interp(
            ggd.electrons.density,
            update_TPS_mats(it, fix_ep_grid_ggd_idx, dd, cell_grid_subset, tps_mats),
            cell_grid_subset,
        )(
            rintersect,
            zintersect,
        )
    end

    return nesep
end
