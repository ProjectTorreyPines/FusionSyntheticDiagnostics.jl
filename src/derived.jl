# File for holding calculations for derived quantities. The inputs should come from
# basic inputs/boundary conditions (like injected power) or from synthetic diagnostic
# outputs (like line average density).

"""
References
[Eldon 2022 PPCF] D. Eldon, et al., Plasma Phys. Control. Fusion 64 (2022) 075002 https://doi.org/10.1088/1361-6587/ac6ff9
[Eich 2013 JNM]
[Eich 2013 NF]
"""

export get_powr_from_dd,
    calc_conducted_loss_power, calc_loss_power, calc_heat_flux_width, calc_q_cyl,
    find_OMP_RZ, read_B_theta_OMP

"""
    function find_OMP_RZ(dd::IMASDD.dd;
        grid_ggd_idx::Int=1,
        space_number::Int=1,
        midplane_grid_subset::Int=11,
        separatrix_grid_subset::Int=16,
    )::Tuple

Reads poloidal magnetic field at the outer midplane separatrix.

  - grid_ggd_idx: index of the grid_ggd set to use within edge_profiles. That is,
    dd.edge_profiles.grid_ggd[grid_ggd_idx] will the source of grid data. To use
    a source other than edge_profiles, use find_OMP_RZ(grid_ggd;) instead.
  - space_number: used for indexing to select different spaces within grid_ggd
  - midplane_grid_subset: index of the subset that contains the 1D outboard midplane
  - separatrix_grid_subset: index of the subset that contains the 1D separatrix
"""
function find_OMP_RZ(dd::IMASDD.dd;
    grid_ggd_idx::Int=1,
    space_number::Int=1,
    midplane_grid_subset::Int=11,
    separatrix_grid_subset::Int=16,
)::Tuple
    grid_ggd = dd.edge_profiles.grid_ggd[grid_ggd_idx]
    return find_OMP_RZ(grid_ggd;
        space_number=space_number,
        midplane_grid_subset=midplane_grid_subset,
        separatrix_grid_subset=separatrix_grid_subset,
    )
end

"""
    function find_OMP_RZ(grid_ggd::IMASDD.edge_profiles__grid_ggd{Float64};
        space_number::Int=1,
        midplane_grid_subset::Int=11,
        separatrix_grid_subset::Int=16,
    )::Tuple

Reads poloidal magnetic field at the outer midplane separatrix.

  - grid_ggd: reference to the grid_ggd item to use for getting geometry data, such as
    dd.equilibrium.grid_ggd[1] or dd.edge_profiles.grid_ggd[1].
  - space_number: used for indexing to select different spaces within grid_ggd
  - midplane_grid_subset: index of the subset that contains the 1D outboard midplane
  - separatrix_grid_subset: index of the subset that contains the 1D separatrix
"""
function find_OMP_RZ(grid_ggd::IMASDD.edge_profiles__grid_ggd{Float64};
    space_number::Int=1,
    midplane_grid_subset::Int=11,
    separatrix_grid_subset::Int=16,
)::Tuple
    space = grid_ggd.space[space_number]
    midplane_sub = get_grid_subset(grid_ggd, midplane_grid_subset)
    separatrix_sub = get_grid_subset(grid_ggd, separatrix_grid_subset)
    intersect_element = subset_do(
        intersect, midplane_sub.element, separatrix_sub.element; space=space,
        use_nodes=true)
    intersect_index = intersect_element[1].object[1].index
    intersect_point = space.objects_per_dimension[1].object[intersect_index].geometry
    R_OMP = intersect_point[1]
    Z_OMP = intersect_point[2]
    return R_OMP, Z_OMP
end

"""
    read_B_theta_OMP(dd::IMASDD.dd;
        grid_ggd_idx::Int=1,
        grid_ggds=nothing,
        space_number::Int=1,
        cell_grid_subset::Int=5,
        midplane_grid_subset::Int=11,
        separatrix_grid_subset::Int=16,
    )::Array{Float64}

Reads poloidal magnetic field at the outer midplane separatrix.

  - dd: an IMAS data dictionary with appropriate data
  - grid_ggd_idx: index of the GGD grid to use. It is often the case that there is only 1
  - grid_ggds: optional reference to the set of grid_ggd items, such as
    dd.equilibrium.grid_ggd or dd.edge_profiles.grid_ggd . This is used in case the
    grid_ggd data are stored separately from equilibrium ggd data as could be the case to
    avoid duplicating the entire grid_ggd.
  - space_number: used for indexing to select different spaces within grid_ggd
  - cell_grid_subset: index of the subset that contains all cells. This should normally
    be 5 for all SOLPS cells, but might be -5 in cases where extended data are added and
    native B2 cells are to be distinguished from extended mesh cells.
  - midplane_grid_subset: index of the subset that contains the 1D outboard midplane
  - separatrix_grid_subset: index of the subset that contains the 1D separatrix
"""
function read_B_theta_OMP(dd::IMASDD.dd;
    grid_ggd_idx::Int=1,
    grid_ggds=nothing,
    space_number::Int=1,
    cell_grid_subset::Int=5,
    midplane_grid_subset::Int=11,
    separatrix_grid_subset::Int=16,
)::Array{Float64}
    # Form empty output array
    nslices = length(dd.equilibrium.time_slice)
    if nslices == 0
        error(
            "No slices in dd.equilibrium.time_slice. ",
            "Cannot read Btheta from empty equilibrium data."
        )
    end
    B_θ_OMP = Array{Float64}(undef, nslices)

    # Get out references to basic ggd objects
    if grid_ggds == nothing
        grid_ggds = dd.equilibrium.grid_ggd
    end
    fix_ep_grid_ggd_idx = length(grid_ggds) == 1
    grid_ggd = grid_ggds[grid_ggd_idx]
    space = grid_ggd.space[space_number]

    rintersect, zintersect = find_OMP_RZ(
        grid_ggd;
        space_number=space_number,
        midplane_grid_subset=midplane_grid_subset,
        separatrix_grid_subset=separatrix_grid_subset,
    )

    tps_mats = get_TPS_mats(grid_ggd, cell_grid_subset)

    for it ∈ 1:nslices
        if length(dd.equilibrium.time_slice[it].ggd) < 1
            B_θ_OMP[it] = read_B_theta_OMP_no_ggd(dd, time_idx=it)
            # error("Equilibrium time slice #$it has no GGD data; cannot proceed.")
        else
            ggd = dd.equilibrium.time_slice[it].ggd[1]

            bz = interp(
                ggd.b_field_z,
                update_TPS_mats(it, fix_ep_grid_ggd_idx, dd, cell_grid_subset, tps_mats),
                cell_grid_subset,
            )(
                rintersect,
                zintersect,
            )
            br = interp(
                ggd.b_field_r,
                update_TPS_mats(it, fix_ep_grid_ggd_idx, dd, cell_grid_subset, tps_mats),
                cell_grid_subset,
            )(
                rintersect,
                zintersect,
            )
            B_θ_OMP[it] = sqrt(bz^2 + br^2) * sign(bz)
        end
    end
    return B_θ_OMP
end

"""
    function read_B_theta_OMP_no_ggd(dd::IMASDD.dd; time_idx::Int=1)

Uses flux map in equilibrium (which is normally on a rectangular grid) to get the
magnetic field at the midplane
"""
function read_B_theta_OMP_no_ggd(dd::IMASDD.dd; time_idx::Int=1)

    # r0 = dd.summary.global_quantities.r0.value
    # b0 = dd.summary.global_quantities.b0.value[time_idx]

    gq = dd.equilibrium.time_slice[time_idx].global_quantities
    psib = gq.psi_boundary
    psia = gq.psi_axis
    zmaxis = gq.magnetic_axis.z

    p2 = dd.equilibrium.time_slice[time_idx].profiles_2d[1]
    psirz = p2.psi
    r = p2.grid.dim1
    z = p2.grid.dim2

    psin = (psirz .- psia) ./ (psib .- psia)
    

    dpsidr, dpsidz = IMASDD.gradient(r, z, psirz)
    br = -dpsidz ./ r
    bz = dpsidr ./ r
    #bphi = r0 * b0 ./ r

    # find closest row to zmaxis
    dz = abs.(z .- zmaxis)
    i1 = argmin(dz)
    dz2 = dz[dz .!= dz[i1]]
    i2 = argmin(abs.(dz .- dz2[argmin(dz2)]))
    wi2 = (z[i1] - zmaxis) / (z[i1] - z[i2])
    wi1 = (zmaxis - z[i2]) / (z[i1] - z[i2])
    psin_omp = psin[:, i1] .* wi1 .+ psin[:, i2] .* wi2
    dp = abs.(psin_omp .- 1)
    j1 = argmin(dp)
    dp2 = dp[dp .!= dp[j1]]
    j2 = argmin(abs.(dp .- dp2[argmin(dp2)]))
    wj2 = (psin_omp[j1] - 1) / (psin_omp[j1] - psin_omp[j2])
    wj1 = (1 - psin_omp[j2]) / (psin_omp[j1] - psin_omp[j2])

    br_omp = (
        br[i1, j1] * wi1 * wj1 + 
        br[i2, j1] * wi2 * wj1 + 
        br[i1, j2] * wi1 * wj2 + 
        br[i2, j2] * wi2 * wj2
    )
    bz_omp = (
        bz[i1, j1] * wi1 * wj1 + 
        bz[i2, j1] * wi2 * wj1 + 
        bz[i1, j2] * wi1 * wj2 + 
        bz[i2, j2] * wi2 * wj2
    )
    
    B_θ_OMP = sqrt(br_omp^2 + bz_omp^2) * sign(bz_omp)
    return B_θ_OMP
end

# """
#     function mag_field_from_equil_to_ggd!(
#         dd::IMASDD.dd; time_idx::Int=1, grid_ggd_source=nothing
#     )

# Uses flux map in equilibrium (which is normally on a rectangular grid) along with a 
# grid_ggd definition instance to define magnetic field components on the ggd mesh.
# """
# function mag_field_from_equil_to_ggd!(
#     dd::IMASDD.dd; time_idx::Int=1, grid_ggd_source=nothing,
# )
#     if grid_ggd_source == nothing
#         grid_ggd_source = dd.edge_profiles.grid_ggd[1]
#         # The other natural source for grid_ggd would be
#         # dd.equilibrium.grid_ggd[1]
#         # But if equilibrium data were already on ggd, we wouldn't need this function at
#         # all. So the most likely use case is anticipated to be using a rectangular flux
#         # map and a GGD from SOLPS via edge_profiles.
#     end

#     r0 = dd.summary.global_quantities.r0.value
#     b0 = dd.summary.global_quantities.b0.value[time_idx]

#     p2 = dd.equilibrium.time_slice[time_idx].profiles_2d
#     psirz = p2.psi
#     r = p2.grid.dim1
#     z = p2.grid.dim2

#     dpsidr, dpsidz = gradient(r, z, psirz)
#     br = -dpsidz ./ r
#     bz = dpsidr ./ r
#     bphi = r0 * b0 ./ r
#
#     ....
# end

"""
    get_power_from_dd(dd::IMASDD.dd)

Utility for extracting power values from the data dictionary. Intended for internal use.
"""
function get_power_from_dd(dd::IMASDD.dd)
    time = dd.summary.time
    W = dd.summary.global_quantities.energy_mhd.value
    P_rad_core = dd.summary.global_quantities.power_radiated_inside_lcfs.value
    P_OHM = dd.summary.global_quantities.power_ohm.value
    P_NBI = dd.summary.heating_current_drive.power_nbi.value
    P_ECH = dd.summary.heating_current_drive.power_ec.value
    P_ICH = dd.summary.heating_current_drive.power_ic.value
    P_LH = dd.summary.heating_current_drive.power_lh.value
    P_fusion = dd.summary.fusion.power.value.value
    P_other = dd.summary.heating_current_drive.power_additional
    return time, W, P_rad_core, P_OHM, P_NBI, P_ECH, P_ICH, P_LH, P_fusion, P_other
end

"""
    set_default_power_arrays(
        nt::Int;
        P_NBI=nothing,
        P_ECH=nothing,
        P_ICH=nothing,
        P_LH=nothing,
        P_fusion=nothing,
        P_other=nothing,
    )

Utility for setting up zero-filled arrays of the same length to cover power values
that weren't provided. Intended for internal use.
"""
function set_default_power_arrays(
    nt::Int;
    P_NBI=nothing,
    P_ECH=nothing,
    P_ICH=nothing,
    P_LH=nothing,
    P_fusion=nothing,
    P_other=nothing,
)
    if P_NBI == nothing
        P_NBI = zeros(nt)
    end
    if P_ECH == nothing
        P_ECH = zeros(nt)
    end
    if P_ICH == nothing
        P_ICH = zeros(nt)
    end
    if P_LH == nothing
        P_LH = zeros(nt)
    end
    if P_fusion == nothing
        P_fusion = zeros(nt)
    end
    if P_other == nothing
        P_other = zeros(nt)
    end
    return P_NBI, P_ECH, P_ICH, P_LH, P_fusion, P_other
end

"""
    calc_conducted_loss_power(
        time::Array{Float64},
        W::Array{Float64},
        P_OHM::Array{Float64};
        P_NBI=nothing,
        P_ECH=nothing,
        P_ICH=nothing,
        P_LH=nothing,
        P_fusion=nothing,
        P_other=nothing,
    )::Array{Float64}

Simple calculation of loss power if only conducted power is considered (neglecting
radiation). This normally would function as an inner layer and be called by the version
of this function that accepts an IMASDD.dd instance.
"""
function calc_conducted_loss_power(
    time::Array{Float64},
    W::Array{Float64},
    P_OHM::Array{Float64};
    P_NBI=nothing,
    P_ECH=nothing,
    P_ICH=nothing,
    P_LH=nothing,
    P_fusion=nothing,
    P_other=nothing,
)::Array{Float64}
    nt = length(time)
    P_NBI, P_ECH, P_ICH, P_LH, P_fusion, P_other = set_default_power_arrays(
        nt;
        P_NBI=P_NBI,
        P_ECH=P_ECH,
        P_ICH=P_ICH,
        P_LH=P_LH,
        P_fusion=P_fusion,
        P_other=P_other,
    )
    dWdt = IMASDD.gradient(time, W)
    P_cond = [
        calc_conducted_loss_power(
            dWdt[i],
            P_OHM[i];
            P_NBI=P_NBI[i],
            P_ECH=P_ECH[i],
            P_ICH=P_ICH[i],
            P_LH=P_LH[i],
            P_fusion=P_fusion[i],
            P_other=P_other[i],
        ) for i ∈ 1:nt
    ]
    return P_cond
end

"""
    function calc_conducted_loss_power(
        dWdt::Float64,
        P_OHM::Float64;
        P_NBI::Float64=0.0,
        P_ECH::Float64=0.0,
        P_ICH::Float64=0.0,
        P_LH::Float64=0.0,
        P_fusion::Float64=0.0,
        P_other::Float64=0.0,
    )::Float64

Simple calculation of loss power if only conducted power is considered (neglecting
radiation). This normally would function as an inner layer and be called by the version
of this function that accepts an IMASDD.dd instance.
"""
function calc_conducted_loss_power(
    dWdt::Float64,
    P_OHM::Float64;
    P_NBI::Float64=0.0,
    P_ECH::Float64=0.0,
    P_ICH::Float64=0.0,
    P_LH::Float64=0.0,
    P_fusion::Float64=0.0,
    P_other::Float64=0.0,
)::Float64
    P_input = P_NBI + P_OHM + P_ECH + P_ICH + P_LH + P_fusion + P_other
    P_cond = P_input - dWdt
    return P_cond
end

"""
    function calc_loss_power(
        time::Array{Float64},
        W::Array{Float64},
        P_rad_core::Array{Float64},
        P_OHM::Array{Float64};
        P_NBI=nothing,
        P_ECH=nothing,
        P_ICH=nothing,
        P_LH=nothing,
        P_fusion=nothing,
        P_other=nothing,
    )::Array{Float64}

Simple calculation of loss power across the separatrix. Core radiation is included. This
normally would function as an inner layer and be called by the version of this function
that accepts an IMASDD.dd instance.
"""
function calc_loss_power(
    time::Array{Float64},
    W::Array{Float64},
    P_rad_core::Array{Float64},
    P_OHM::Array{Float64};
    P_NBI=nothing,
    P_ECH=nothing,
    P_ICH=nothing,
    P_LH=nothing,
    P_fusion=nothing,
    P_other=nothing,
)::Array{Float64}
    dWdt = IMASDD.gradient(time, W)
    nt = length(time)
    P_NBI, P_ECH, P_ICH, P_LH, P_fusion, P_other = set_default_power_arrays(
        nt;
        P_NBI=P_NBI,
        P_ECH=P_ECH,
        P_ICH=P_ICH,
        P_LH=P_LH,
        P_fusion=P_fusion,
        P_other=P_other,
    )
    P_cond = [
        calc_conducted_loss_power(
            dWdt[i],
            P_OHM[i];
            P_NBI=P_NBI[i],
            P_ECH=P_ECH[i],
            P_ICH=P_ICH[i],
            P_LH=P_LH[i],
            P_fusion=P_fusion[i],
            P_other=P_other[i],
        ) for i ∈ 1:nt
    ]
    P_SOL = P_cond - P_rad_core
    return P_SOL
end

"""
    function calc_loss_power(
        dWdt::Float64,
        P_rad_core::Float64,
        P_OHM::Float64;
        P_NBI::Float64=0.0,
        P_ECH::Float64=0.0,
        P_ICH::Float64=0.0,
        P_LH::Float64=0.0,
        P_fusion::Float64=0.0,
        P_other::Float64=0.0,
    )::Float64

Simple calculation of loss power across the separatrix. Core radiation is included. This
normally would function as an inner layer and be called by the version of this function
that accepts an IMASDD.dd instance.
"""
function calc_loss_power(
    dWdt::Float64,
    P_rad_core::Float64,
    P_OHM::Float64;
    P_NBI::Float64=0.0,
    P_ECH::Float64=0.0,
    P_ICH::Float64=0.0,
    P_LH::Float64=0.0,
    P_fusion::Float64=0.0,
    P_other::Float64=0.0,
)::Float64
    P_cond = calc_conducted_loss_power(
        dWdt,
        P_OHM;
        P_NBI=P_NBI,
        P_ECH=P_ECH,
        P_ICH=P_ICH,
        P_LH=P_LH,
        P_fusion=P_fusion,
        P_other=P_other,
    )
    P_SOL = P_cond - P_rad_core
    return P_SOL
end

"""
    function calc_loss_power(dd::IMASDD.dd)::Array{Float64}

Simple calculation of loss power across the separatrix. Core radiation is included.
"""
function calc_loss_power(dd::IMASDD.dd)::Array{Float64}
    time, W, P_rad_core, P_OHM, P_NBI, P_ECH, P_ICH, P_LH, P_fusion, P_other =
        get_power_from_dd(dd)
    P_SOL = calc_loss_power(
        time,
        W,
        P_rad_core,
        P_OHM;
        P_NBI=P_NBI,
        P_ECH=P_ECH,
        P_ICH=P_ICH,
        P_LH=P_LH,
        P_fusion=P_fusion,
        P_other=P_other,
    )
    return P_SOL
end

"""
    calc_q_cyl(
        B_ϕ_axis::Float64,
        Iₚ::Float64,
        aₘᵢₙₒᵣ::Float64,
        R_geo::Float64,
        κ::Float64,
    )::Float64

Calculates the cylinderical version of safety factor.
B_ϕ_axis: the toroidal magnetic field component at the magnetic axis / T
Iₚ: the total plasma current / A
aₘᵢₙₒᵣ: the minor radius of the plasma / m
R_geo: the geometric major radius of the plasma
(the average of R at the two points where the separatrix crosses the midplane) / m
κ: elongation / unitless
"""
function calc_q_cyl(
    B_ϕ_axis::Float64,
    Iₚ::Float64,
    aₘᵢₙₒᵣ::Float64,
    R_geo::Float64,
    κ::Float64,
)::Float64
    # Copied from Equation 19 of [Eldon 2022 PPCF] which is from Equation 6 of [Eich 2013 JNM]
    ε = aₘᵢₙₒᵣ / R_geo
    μ₀ = π * 4e-7  # H / m
    q_cyl = π * aₘᵢₙₒᵣ * ε * B_ϕ_axis / (μ₀ * Iₚ) * (1 + κ^2)
    return q_cyl  # unitless
end

"""
    calc_heat_flux_width(dd::IMASDD.dd; version::Int=1)

Calculates heat flux width from a scaling law. Different regressions are available
by selecting different versions of the scaling law.

Version 1 takes the form of equation 5 of [Eich 2013 JNM] with values from the
top row of table 6 from [Eich 2013 NF]
Version 2 is a simple constant times the poloidal field at the outboard midplane
from regression #14 in table 3 of [Eich 2013 NF]
"""
function calc_heat_flux_width(dd::IMASDD.dd; version::Int=1)::Array{Float64}
    if version == 1
        B_ϕ_axis =
            dd.equilibrium.time_slice[:].global_quantities.magnetic_axis.b_field_tor
        P_SOL = calc_loss_power(dd)
        R_geo = dd.summary.boundary.geometric_axis_r.value
        aₘᵢₙₒᵣ = dd.summary.boundary.minor_radius.value
        κ = dd.summary.boundary.elongation.value
        Iₚ = dd.summary.global_quantities.ip.value
        λq = calc_heat_flux_width.(B_ϕ_axis, P_SOL, R_geo, aₘᵢₙₒᵣ, κ, Iₚ)  # mm
    elseif version == 2
        B_θ_OMP =
            λq = calc_heat_flux_width.(B_θ_OMP)  # mm
    end
    return λq  # mm
end

"""
    function calc_heat_flux_width(
        B_ϕ_axis::Float64,
        P_SOL::Float64,
        R_geo::Float64,
        aₘᵢₙₒᵣ::Float64,
        κ::Float64,
        Iₚ::Float64,
    )::Float64

Heat flux width from scaling law, using the form of equation 5 of [Eich 2013 JNM] with
values from the top row of table 6 from [Eich 2013 NF]. This one accepts more base
quantities and calculates q_cyl.

B_ϕ_axis: the toroidal magnetic field component at the magnetic axis / T
P_SOL: Power crossing the last closed flux surface / MW
R_geo: the geometric major radius of the plasma
(the average of R at the two points where the separatrix crosses the midplane) / m
aₘᵢₙₒᵣ: the minor radius of the plasma / m
κ: elongation / unitless
Iₚ: the total plasma current / A
"""
function calc_heat_flux_width(
    B_ϕ_axis::Float64,
    P_SOL::Float64,
    R_geo::Float64,
    aₘᵢₙₒᵣ::Float64,
    κ::Float64,
    Iₚ::Float64,
)::Float64
    q_cyl = calc_q_cyl(B_ϕ_axis, Iₚ, aₘᵢₙₒᵣ, R_geo, κ)  # Unitless
    λq = calc_heat_flux_width(B_ϕ_axis, P_SOL, R_geo, q_cyl)  # mm
    return λq
end

"""
    function calc_heat_flux_width(
        B_ϕ_axis::Float64,
        P_SOL::Float64,
        R_geo::Float64,
        q_cyl::Float64,
    )::Float64

Heat flux width from scaling law, using the form of equation 5 of [Eich 2013 JNM] with
values from the top row of table 6 from [Eich 2013 NF]. This one needs q_cyl to be
evaluated first.

B_ϕ_axis: the toroidal magnetic field component at the magnetic axis / T
P_SOL: Power crossing the last closed flux surface / MW
R_geo: the geometric major radius of the plasma
(the average of R at the two points where the separatrix crosses the midplane) / m
q_cyl: cylindrical safety factor / unitless
"""
function calc_heat_flux_width(
    B_ϕ_axis::Float64,
    P_SOL::Float64,
    R_geo::Float64,
    q_cyl::Float64,
)::Float64
    # Copied from Equation 18 of [Eldon 2022 PPCF] which came from [Eich 2013 JNM] and [Eich 2013 NF].
    # Specifically, values from table 6 (top row for JET/DIII-D/AUT) of [Eich 2013 NF] were used with
    # equation 5 of [Eich 2013 JNM]
    λq =
        0.86 * abs(B_ϕ_axis)^(-0.8) * abs(q_cyl)^(1.11) * abs(P_SOL)^0.11 *
        abs(R_geo)^(-0.13)  # mm
    return λq
end

"""
    calc_heat_flux_width(B_θ_OMP::Float64)::Float64

Heat flux width from scaling law, using a simple constant times the poloidal field at the
outboard midplane from regression #14 in table 3 of [Eich 2013 NF]

B_θ_OMP: poloidal magnetic field at the outboard midplane separatrix
"""
function calc_heat_flux_width(B_θ_OMP::Float64)::Float64
    # From Table 3 of [Eich 2013 NF], row 3, regression #14 for all tokamaks
    λq = 0.63 * abs(B_θ_OMP)^(-1.19)  # mm
    return λq
end
