import PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e
import SD4SOLPS: fill_in_extrapolated_core_profile!, check_rho_1d,
    add_rho_to_equilibrium!, core_profile_2d
import QuadGK: quadgk, BatchIntegrand
using GGDUtils
using SOLPS2IMAS: SOLPS2IMAS
using ProgressBars: ProgressBar

default_ifo = "$(@__DIR__)/default_interferometer.json"

"""
    add_interferometer(
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
    config::String=default_ifo,

)::OMAS.dd

Add interferometer to IMAS structure using a JSON file and compute the
line integrated electron density if not present
"""
function add_interferometer!(config::String=default_ifo,
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
)::OMAS.dd
    if endswith(config, ".json")
        OMAS.json2imas(config, ids)
    else
        error("Only JSON files are supported.")
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
function add_interferometer!(config::Dict{Symbol, Any},
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
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
function compute_interferometer(@nospecialize(ids::OMAS.dd), rtol::Float64=1e-1)
    fix_eq_time_idx = length(ids.equilibrium.time_slice) == 1
    fix_ep_grid_ggd_idx = length(ids.edge_profiles.grid_ggd) == 1

    ep_grid_ggd = ids.edge_profiles.grid_ggd[1]
    ep_space = ep_grid_ggd.space[1]
    all_cells = SOLPS2IMAS.get_grid_subset_with_index(ep_grid_ggd, 5)
    SOLPS_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
    SOLPS_bnd.element = SOLPS2IMAS.get_subset_boundary(ep_space, all_cells)
    core_bnd = SOLPS2IMAS.get_grid_subset_with_index(ep_grid_ggd, 15)
    core = SOLPS2IMAS.get_grid_subset_with_index(ep_grid_ggd, 22)
    sol = SOLPS2IMAS.get_grid_subset_with_index(ep_grid_ggd, 23)

    cp_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
    cp_bnd.element =
        SOLPS2IMAS.subset_do(
            intersect,
            SOLPS2IMAS.get_subset_boundary(ep_space, sol),
            SOLPS2IMAS.get_subset_boundary(ep_space, core),
        )
    SOLPS_out_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
    SOLPS_out_bnd.element =
        SOLPS2IMAS.subset_do(setdiff, SOLPS_bnd.element, core_bnd.element)
    ep_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
    ep_bnd.element =
        SOLPS2IMAS.subset_do(union, SOLPS_out_bnd.element, cp_bnd.element)

    TPS_mats_all_cells = get_TPS_mats(ep_grid_ggd, 5)

    ep_n_e_list = []
    cp_n_e_list = []

    for ch ∈ ids.interferometer.channel
        k = [2 * π / ch.wavelength[ii].value for ii ∈ 1:2]
        for i1 ∈ 1:2
            lam = ch.wavelength[i1]
            i2 = i1 % 2 + 1
            if lam.phase_to_n_e_line == 0
                # Taken from https://doi.org/10.1063/1.1138037
                lam.phase_to_n_e_line =
                    (
                        2 * m_e * ε_0 * c_0^2 / e^2 * (k[i1] * k[i2]^2) /
                        (k[i2]^2 - k[i1]^2)
                    ).val
            end
        end
        # Special case when measurement has not been made but edge profile data exists
        if length(ch.n_e_line.time) == 0 && length(ids.edge_profiles.ggd) > 0
            epggd = ids.edge_profiles.ggd
            cpp1d = ids.core_profiles.profiles_1d
            nt = length(epggd)
            ch.n_e_line.time = zeros(nt)
            ch.n_e_line_average.time = zeros(nt)
            ch.n_e_line.data = zeros(nt)
            ch.n_e_line_average.data = zeros(nt)
            for lam ∈ ch.wavelength
                lam.phase_corrected.time = zeros(nt)
                lam.phase_corrected.data = zeros(nt)
            end
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
            # Parametrize line of sight
            fp = rzphi2xyz(ch.line_of_sight.first_point)
            sp = rzphi2xyz(ch.line_of_sight.second_point)
            tp = rzphi2xyz(ch.line_of_sight.third_point)
            if ch.line_of_sight.third_point ==
               OMAS.interferometer__channel___line_of_sight__third_point()
                los = (s::Float64) -> xyz2rz((fp .+ s .* (sp .- fp))...)
                dl_ds = (s::Float64) -> sqrt(sum((sp .- fp) .^ 2))
                chord_points = (fp, sp)
            else
                los =
                    (s::Float64) ->
                        if s <= 0.5
                            return xyz2rz((fp .+ 2 * s .* (sp .- fp))...)
                        else
                            return xyz2rz((sp .+ 2 * (s - 0.5) .* (tp .- sp))...)
                        end
                dl_ds = (s::Float64) -> if s <= 0.5
                    return sqrt(sum((2 .* (sp .- fp)) .^ 2))
                else
                    return sqrt(sum((2 .* (sp .- tp)) .^ 2))
                end
                chord_points = (fp, sp, tp)
            end

            # ep_segs = get_intersections(ep_bnd, ep_space, chord_points, los)
            # cp_segs = get_intersections(cp_bnd, ep_space, chord_points, los)

            for ii ∈ eachindex(epggd)
                ch.n_e_line.time[ii] = epggd[ii].time
                ch.n_e_line_average.time[ii] = epggd[ii].time
                eq_time_idx = fix_eq_time_idx ? 1 : ii
                if !fix_ep_grid_ggd_idx
                    ep_grid_ggd = ids.edge_profiles.grid_ggd[ii]
                    ep_space = ep_grid_ggd.space[1]
                    all_cells = SOLPS2IMAS.get_grid_subset_with_index(ep_grid_ggd, 5)
                    SOLPS_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
                    SOLPS_bnd.element =
                        SOLPS2IMAS.get_subset_boundary(ep_space, all_cells)
                    core_bnd = SOLPS2IMAS.get_grid_subset_with_index(ep_grid_ggd, 15)
                    core = SOLPS2IMAS.get_grid_subset_with_index(ep_grid_ggd, 22)
                    sol = SOLPS2IMAS.get_grid_subset_with_index(ep_grid_ggd, 23)

                    cp_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
                    cp_bnd.element =
                        SOLPS2IMAS.subset_do(
                            intersect,
                            SOLPS2IMAS.get_subset_boundary(ep_space, sol),
                            SOLPS2IMAS.get_subset_boundary(ep_space, core),
                        )
                    SOLPS_out_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
                    SOLPS_out_bnd.element =
                        SOLPS2IMAS.subset_do(
                            setdiff,
                            SOLPS_bnd.element,
                            core_bnd.element,
                        )
                    ep_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
                    ep_bnd.element =
                        SOLPS2IMAS.subset_do(
                            union,
                            SOLPS_out_bnd.element,
                            cp_bnd.element,
                        )

                    TPS_mats_all_cells = get_TPS_mats(ep_grid_ggd, 5)

                    ep_segs = get_intersections(ep_bnd, ep_space, chord_points, los)
                    cp_segs = get_intersections(cp_bnd, ep_space, chord_points, los)
                end

                if length(ep_n_e_list) < ii
                    append!(
                        ep_n_e_list,
                        [interp(epggd[ii].electrons.density, TPS_mats_all_cells, 5)],
                    )
                end
                if length(cp_n_e_list) < ii
                    cp_p1d = ids.core_profiles.profiles_1d[ii]
                    eqt = ids.equilibrium.time_slice[eq_time_idx]
                    append!(
                        cp_n_e_list,
                        [interp(cp_p1d.electrons.density, cp_p1d, eqt)],
                    )
                end

                ep_n_e = ep_n_e_list[ii]
                cp_n_e = cp_n_e_list[ii]

                # integ_ep = (s) -> ep_n_e(los(s)...) * dl_ds(s)
                # integ_cp = (s) -> cp_n_e.(los(s)...) * dl_ds(s)

                # int_n_e = 0
                # # println("In Edge Profiles:")
                # for seg ∈ ep_segs
                #     # println(seg[1], " to ", seg[2])
                #     # println("Slope: ", dl_ds((seg[1] + seg[2]) / 2))
                #     # println("Inegrand: ", integ_ep((seg[1] + seg[2]) / 2))
                #     seg_int = quadgk(integ_ep, seg[1], seg[2]; rtol=rtol)[1]
                #     int_n_e += seg_int
                #     # println("This segment: ", seg_int, "Total: ", int_n_e)
                # end
                # # println("In Core Profiles:")
                # for seg ∈ cp_segs
                #     # println(seg[1], " to ", seg[2])
                #     # println("Slope: ", dl_ds((seg[1] + seg[2]) / 2))
                #     # println("Inegrand: ", integ_cp((seg[1] + seg[2]) / 2))
                #     seg_int = quadgk(integ_cp, seg[1], seg[2]; rtol=rtol)[1]
                #     int_n_e += seg_int
                #     # println("This segment: ", seg_int, "Total: ", int_n_e)
                # end
                # ch.n_e_line.data[ii] = int_n_e
                integrand(s) =
                    n_e(
                        cp_n_e,
                        ep_n_e,
                        core_bnd,
                        SOLPS_bnd,
                        ep_space,
                        los(s)...,
                    ) .* dl_ds(s)
                # function integrand!(y, s)
                #     n = Threads.nthreads()
                #     Threads.@threads for i ∈ 1:n
                #         y[i:n:end] .= integrand.(@view(s[i:n:end]))
                #     end
                # end
                # println("Calculating integrated electron density")
                # ch.n_e_line.data[ii] = 0.0
                # for seg ∈ valid_segs
                #     @time ch.n_e_line.data[ii] +=
                #         quadgk(
                #             BatchIntegrand{Float64}(integrand!),
                #             seg[1],
                #             seg[2],
                #             rtol,
                #         )[1]
                # end
                ch.n_e_line.data[ii] = quadgk(integrand, 0, 1; rtol=rtol)[1]
                # println(ch.n_e_line.data[ii])
                chord_in_core = get_intersections(cp_bnd, ep_space, chord_points, los)
                core_chord_length = 0.0
                for seg ∈ chord_in_core
                    core_chord_length +=
                        (seg[2] - seg[1]) * dl_ds((seg[1] + seg[2]) / 2)
                end
                # core_chord_integrand(s) =
                #     is_in_core( cp_bnd, ep_space, los(s)...) .* dl_ds(s)
                # println("Calculating core chord length")
                # @time core_chord_length = quadgk(core_chord_integrand, 0, 1, rtol)[1]
                # println(core_chord_length)
                ch.n_e_line_average.data[ii] =
                    ch.n_e_line.data[ii] / core_chord_length
                # println("Average: ", ch.n_e_line_average.data[ii])
                for lam ∈ ch.wavelength
                    lam.phase_corrected.time[ii] = epggd[ii].time
                    lam.phase_corrected.data[ii] =
                        ch.n_e_line.data[ii] / lam.phase_to_n_e_line
                end
            end
        end
    end
end

function rzphi2xyz(
    point::Union{OMAS.interferometer__channel___line_of_sight__first_point,
        OMAS.interferometer__channel___line_of_sight__second_point,
        OMAS.interferometer__channel___line_of_sight__third_point},
)
    r, z, phi = point.r, point.z, point.phi
    return r * cos(phi), r * sin(phi), z
end

function xyz2rz(x::Float64, y::Float64, z::Float64)
    r = sqrt(x^2 + y^2)
    return r, z
end

function n_e(
    cp_n_e::Function,
    ep_n_e::Function,
    core_bnd::OMAS.edge_profiles__grid_ggd___grid_subset,
    SOLPS_bnd::OMAS.edge_profiles__grid_ggd___grid_subset,
    ep_space::OMAS.edge_profiles__grid_ggd___space,
    r::Float64,
    z::Float64,
)
    if (r, z) ∈ (core_bnd, ep_space)
        return cp_n_e(r, z)
    else# if (r, z) ∈ (SOLPS_bnd, ep_space)
        return ep_n_e(r, z)
        # else
        #     return 0.0
    end
end

function get_intersections(subset, space, points, los)
    nodes = space.objects_per_dimension[1].object
    edges = space.objects_per_dimension[2].object
    if length(points) == 2
        fp, sp = points
        three_points = false
    else
        fp, sp, tp = points
        three_points = true
    end
    segs = [0.0]
    for ele ∈ subset.element
        edge = edges[ele.object[1].index]
        e1 = Tuple(nodes[edge.nodes[1]].geometry)
        e2 = Tuple(nodes[edge.nodes[2]].geometry)
        s, s2 = intersection_s(e1, e2, fp, sp)
        if 0 <= s < 1 && 0 <= s2 < 1
            if three_points
                append!(segs, s / 2)
                s, s2 = intersection_s(e1, e2, sp, tp)
                if 0 <= s < 1 && 0 <= s2 < 1
                    append!(segs, s / 2 + 0.5)
                end
            else
                append!(segs, s)
            end
        end
    end
    append!(segs, 1.0)
    sort!(segs)
    filt_segs = Tuple{Float64, Float64}[]
    for ii ∈ 1:length(segs)-1
        check_at = (segs[ii] + segs[ii+1]) / 2
        # print("Seg : ", los(segs[ii]), " to ", los(segs[ii+1]), ", Mid: ")
        # print(check_at, ": ", los(check_at), ": ")
        if los(check_at) ∈ (subset, space)
            append!(filt_segs, [(segs[ii], segs[ii+1])])
            #     println("Present")
            # else
            #     println("Absent")
        end
    end
    return filt_segs
end

function intersection_s(
    fp::Tuple{Float64, Float64},
    sp::Tuple{Float64, Float64},
    l1::Tuple{Float64, Float64},
    l2::Tuple{Float64, Float64},
)
    den1 = (sp[1] - fp[1]) * (l2[2] - l1[2]) - (sp[2] - fp[2]) * (l2[1] - l1[1])
    num1 = (sp[1] - fp[1]) * (fp[2] - l1[2]) - (sp[2] - fp[2]) * (fp[1] - l1[1])
    den2 = (l2[1] - l1[1]) * (sp[2] - fp[2]) - (l2[2] - l1[2]) * (sp[1] - fp[1])
    num2 = (l2[1] - l1[1]) * (l1[2] - fp[2]) - (l2[2] - l1[2]) * (l1[1] - fp[1])
    return num1 / den1, num2 / den2
end

function intersection_s(
    fp::Tuple{Float64, Float64, Float64},
    sp::Tuple{Float64, Float64, Float64},
    l1::Tuple{Float64, Float64},
    l2::Tuple{Float64, Float64},
)
    return intersection_s(xyz2rz(fp...), xyz2rz(sp...), l1, l2)
end

function intersection_s(
    fp::Tuple{Float64, Float64},
    sp::Tuple{Float64, Float64},
    l1::Tuple{Float64, Float64, Float64},
    l2::Tuple{Float64, Float64, Float64},
)
    return intersection_s(fp, sp, xyz2rz(l1...), xyz2rz(l2...))
end
