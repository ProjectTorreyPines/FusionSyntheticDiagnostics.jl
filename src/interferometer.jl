import PhysicalConstants.CODATA2018: c_0, ε_0, m_e, e
import QuadGK: quadgk, BatchIntegrand
import GGDUtils:
    interp, get_grid_subset_with_index, get_subset_boundary, subset_do, get_TPS_mats

default_ifo = "$(@__DIR__)/default_interferometer.json"

"""
    add_interferometer(
    @nospecialize(ids::OMAS.dd)=OMAS.dd(),
    config::String=default_ifo,

)::OMAS.dd

Add interferometer to IMAS structure using a JSON file and compute the
line integrated electron density if not present
"""
function add_interferometer!(
    config::String=default_ifo,
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
function add_interferometer!(
    config::Dict{Symbol, Any},
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
function compute_interferometer(@nospecialize(ids::OMAS.dd), rtol::Float64=1e-3)
    fix_eq_time_idx = length(ids.equilibrium.time_slice) == 1
    fix_ep_grid_ggd_idx = length(ids.edge_profiles.grid_ggd) == 1

    ep_grid_ggd = ids.edge_profiles.grid_ggd[1]
    ep_space = ep_grid_ggd.space[1]
    sep_bnd, core_bnd = get_sep_core_bnd(ep_grid_ggd)
    # Using -5 for now to use the SOLPS edge profile grid only
    TPS_mats_all_cells = get_TPS_mats(ep_grid_ggd, -5)

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

            # Check if core_profile is available
            if length(cpp1d) != length(epggd)
                error(
                    "Number of edge profiles time slices does not match number of " \
                    "core profile time slices. Please ensure core profile data for " \
                    "electron density is present in the data structure. You might " \
                    "want to run SD4SOLPS.fill_in_extrapolated_core_profile!" \
                    "(dd, \"electrons.density\"; method=core_method))",
                )
            end
            # Parametrize line of sight
            fp = rzphi2xyz(ch.line_of_sight.first_point)
            sp = rzphi2xyz(ch.line_of_sight.second_point)
            tp = rzphi2xyz(ch.line_of_sight.third_point)
            if ch.line_of_sight.third_point ==
               OMAS.interferometer__channel___line_of_sight__third_point()
                chord_points = (fp, sp)
            else
                chord_points = (fp, sp, tp)
            end
            los, dl_ds = get_line_of_sight(chord_points)
            core_chord_length = get_core_chord_length(sep_bnd, ep_space, chord_points)

            for ii ∈ eachindex(epggd)
                ch.n_e_line.time[ii] = epggd[ii].time
                ch.n_e_line_average.time[ii] = epggd[ii].time
                eq_time_idx = fix_eq_time_idx ? 1 : ii
                if !fix_ep_grid_ggd_idx
                    # If grid_ggd is evolving with time, update boundaries
                    ep_grid_ggd = ids.edge_profiles.grid_ggd[ii]
                    ep_space = ep_grid_ggd.space[1]
                    sep_bnd, core_bnd = get_sep_core_bnd(ep_grid_ggd)
                    TPS_mats_all_cells = get_TPS_mats(ep_grid_ggd, 5)
                    core_chord_length =
                        get_core_chord_length(sep_bnd, ep_space, chord_points)
                end

                # Note the grid_subset index 5 is used here to get the electron density
                # data from SOLPS edge profile. This is a bug in SD4SOLPS. On adding
                # edge extension, it should update all quantitites that refered to
                # grid_subset index 5 to -5.
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

                function integrand(s::Real)
                    r, z = los(s)
                    if (r, z) ∈ (sep_bnd, ep_space)
                        return cp_n_e(r, z) .* dl_ds(s)
                    else
                        return ep_n_e(r, z) .* dl_ds(s)
                    end
                end

                ch.n_e_line.data[ii] = quadgk(integrand, 0, 1; rtol=rtol)[1]
                ch.n_e_line_average.data[ii] = ch.n_e_line.data[ii] / core_chord_length
                for lam ∈ ch.wavelength
                    lam.phase_corrected.time[ii] = epggd[ii].time
                    lam.phase_corrected.data[ii] =
                        ch.n_e_line.data[ii] / lam.phase_to_n_e_line
                end
            end
        end
    end
end

function get_sep_core_bnd(ep_grid_ggd)
    ep_space = ep_grid_ggd.space[1]
    core_bnd = get_grid_subset_with_index(ep_grid_ggd, 15)
    core = get_grid_subset_with_index(ep_grid_ggd, 22)
    sol = get_grid_subset_with_index(ep_grid_ggd, 23)
    sep_bnd = OMAS.edge_profiles__grid_ggd___grid_subset()
    sep_bnd.element =
        subset_do(
            intersect,
            get_subset_boundary(ep_space, sol),
            get_subset_boundary(ep_space, core),
        )
    return sep_bnd, core_bnd
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

function get_line_of_sight(points)
    if length(points) == 2
        fp, sp = points
        los = (s::Float64) -> xyz2rz((fp .+ s .* (sp .- fp))...)
        dl_ds = (s::Float64) -> sqrt(sum((sp .- fp) .^ 2))
    else
        fp, sp, tp = points
        los =
            (s::Float64) ->
                if s <= 0.5
                    return xyz2rz((fp .+ 2 * s .* (sp .- fp))...)
                else
                    return xyz2rz((sp .+ 2 * (s - 0.5) .* (tp .- sp))...)
                end
        dl_ds =
            (s::Float64) ->
                if s <= 0.5
                    return sqrt(sum((2 .* (sp .- fp)) .^ 2))
                else
                    return sqrt(sum((2 .* (sp .- tp)) .^ 2))
                end
    end
    return los, dl_ds
end

function get_intersections(subset, space, points)
    los, dl_ds = get_line_of_sight(points)
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
        if los(check_at) ∈ (subset, space)
            append!(filt_segs, [(segs[ii], segs[ii+1])])
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

function get_core_chord_length(sep_bnd, ep_space, chord_points)
    los, dl_ds = get_line_of_sight(chord_points)
    chord_in_core = get_intersections(sep_bnd, ep_space, chord_points)
    core_chord_length = 0.0
    for seg ∈ chord_in_core
        center_of_seg = (seg[1] + seg[2]) / 2
        core_chord_length += (seg[2] - seg[1]) * dl_ds(center_of_seg)
    end
    return core_chord_length
end
