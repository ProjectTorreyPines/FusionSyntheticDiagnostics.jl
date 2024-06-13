import PhysicalConstants.CODATA2018: c_0, ε_0, m_e, m_u, e
import QuadGK: quadgk, BatchIntegrand
import GGDUtils: interp, get_grid_subset, get_subset_boundary, subset_do, get_TPS_mats
using StaticArrays
using LinearAlgebra
using CoordinateTransformations
using Rotations

export add_bolometer!, compute_bolometer!

"""
    FoV

A structure to represent the field of view of a bolometer channel. It holds the half
angle of the conical field of view and the transformation to go from the field of view
frame to the R, Z, Phi frame. In the field of view frame, the z-axis is along the
conical axis and the origin is the apex of the cone. The positive z direction is away
from the detector and towards the plasma.

Constructors:

    FoV(ha::Float64, fov2RZP::AffineMap{T, U})
     where {T <: Rotation{3, Float64}, U <: SVector{3, Float64}}
"""
mutable struct FoV
    var"ha"::Float64 # Half angle of conical field of view in radians
    var"fov2RZP"::AffineMap{
        T,
        U,
    } where {T <: Rotation{3, Float64}, U <: SVector{3, Float64}}
end

default_bolometer = "$(@__DIR__)/default_bolometer.json"

"""
    add_bolometer!(
        config::Union{String, Dict{Symbol, Any}}=default_bolometer,
        @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
        overwrite::Bool=false, verbose::Bool=false, kwargs...,
    )::IMASDD.dd

Add bolometer to IMAS structure using a `JSON` file or Julia `Dict` and compute the
bolometer outputs. `kwargs` are passed to [`compute_bolometer!`](@ref).
"""
function add_bolometer!(
    config::Union{String, Dict{Symbol, Any}}=default_bolometer,
    @nospecialize(ids::IMASDD.dd)=IMASDD.dd();
    overwrite::Bool=false, verbose::Bool=false, kwargs...,
)::IMASDD.dd
    add_diagnostic!(config, :bolometer, ids; overwrite=overwrite, verbose=verbose)
    compute_bolometer!(ids; kwargs...)
    return ids
end

"""
    compute_bolometer!(
        @nospecialize(ids::IMASDD.dd);
        rtol::Float64=1e-3,
        n_e_gsi::Int=5,
    )

Computed the line integrated electron density from the bolometer data present in
IDS structure for all the chords. The computation is based on the edge profile data
and core profile data present in the IDS structure.
"""
function compute_bolometer!(
    @nospecialize(ids::IMASDD.dd);
    c2o_nop::Int64=12,
)
    FoVs = Dict{String, FoV}()
    for ch ∈ ids.bolometer.channel
        # List of transformation to go from aperture frame to R, Z, Phi frame
        ap2RZP = Array{LinearMap{RotMatrix3{Float64}}}(undef, length(ch.aperture))

        # Add outline to apertures and detector if not present
        for (ii, ap) ∈ enumerate(ch.aperture)
            if ap.geometry_type != 1
                if ap.geometry_type == 2
                    resize!(ap.outline.x1, circle_to_outline_nop)
                    resize!(ap.outline.x2, circle_to_outline_nop)
                    ap.outline.x1 = ap.radius .* cos.(range(0, 2π, c2o_nop))
                    ap.outline.x2 = ap.radius .* sin.(range(0, 2π, c2o_nop))
                    if IMASDD.ismissing(ap, :surface)
                        ap.surface = π * ap.radius^2
                    end
                elseif ap.geometry_type == 3
                    resize!(ap.outline.x1, 4)
                    resize!(ap.outline.x2, 4)
                    ap.outline.x1 = [
                        -ap.x1_width / 2,
                        -ap.x1_width / 2,
                        ap.x1_width / 2,
                        ap.x1_width / 2,
                    ]
                    ap.outline.x2 = [
                        -ap.x2_width / 2,
                        ap.x2_width / 2,
                        ap.x2_width / 2,
                        -ap.x2_width / 2,
                    ]
                    if IMASDD.ismissing(ap, :surface)
                        ap.surface = ap.width * ap.height
                    end
                end
            end
            x1 = SVector{3, Float64}(
                ap.x1_unit_vector.x,
                ap.x1_unit_vector.y,
                ap.x1_unit_vector.z,
            )
            x2 = SVector{3, Float64}(
                ap.x2_unit_vector.x,
                ap.x2_unit_vector.y,
                ap.x2_unit_vector.z,
            )
            x3 = SVector{3, Float64}(
                ap.x3_unit_vector.x,
                ap.x3_unit_vector.y,
                ap.x3_unit_vector.z,
            )
            ap2RZP[ii] = inv(LinearMap(RotMatrix3([x1; x2; x3])))
        end

        det = ch.detector
        if det.geometry_type != 1
            if det.geometry_type == 2
                resize!(det.outline.x1, circle_to_outline_nop)
                resize!(det.outline.x2, circle_to_outline_nop)
                det.outline.x1 = det.radius .* cos.(range(0, 2π, c2o_nop))
                det.outline.x2 = det.radius .* sin.(range(0, 2π, c2o_nop))
                if IMASDD.ismissing(det, :surface)
                    det.surface = π * det.radius^2
                end
            elseif det.geometry_type == 3
                resize!(det.outline.x1, 4)
                resize!(det.outline.x2, 4)
                det.outline.x1 = [
                    -det.x1_width / 2,
                    -det.x1_width / 2,
                    det.x1_width / 2,
                    det.x1_width / 2,
                ]
                det.outline.x2 = [
                    -det.x2_width / 2,
                    det.x2_width / 2,
                    det.x2_width / 2,
                    -det.x2_width / 2,
                ]
                if IMASDD.ismissing(det, :surface)
                    det.surface = det.width * det.height
                end
            end
        end
        x1 = SVector{3, Float64}(
            det.x1_unit_vector.x,
            det.x1_unit_vector.y,
            det.x1_unit_vector.z,
        )
        x2 = SVector{3, Float64}(
            det.x2_unit_vector.x,
            det.x2_unit_vector.y,
            det.x2_unit_vector.z,
        )
        x3 = SVector{3, Float64}(
            det.x3_unit_vector.x,
            det.x3_unit_vector.y,
            det.x3_unit_vector.z,
        )
        # Transformation to go from detector frame to R, Z, Phi frame
        det2RZP = inv(LinearMap(RotMatrix3([x1; x2; x3])))

        # Transformation to go from aperture frame to detector frame
        R_ap2det = Array{LinearMap{RotMatrix3{Float64}}}(undef, length(ch.aperture))
        for (ii, R) ∈ enumerate(R_ap)
            R_ap2det[ii] = inv(det2RZP) ∘ R
        end

        FoVs[ch.identifier] = get_FoV(ch, ap2RZP[end], det2RZP)
    end
end

"""
    get_FoV(
        ch::IMASDD.bolometer_channel,
        ap2RZP::LinearMap{RotMatrix3{Float64}},
        det2RZP::LinearMap{RotMatrix3{Float64}},
    )::FoV

Assuming aperture and detector are arrays of points described in 3D space on same
coordinate axes and that the array of points for aperture lie in a plane and the array
of points for the detector lie in a plane. This function computes the field of view
of the detector from the aperture and returns a tuple of the half angle of the field of
view and the transformation that maps from field of view coordinates to the input
coordinate frame. In the field of view coordinate frame, the field of view is a cone
with the vertex at the origin and the axis along the z-axis, and half angle that is
returned.
"""
function get_FoV(
    ch::IMASDD.bolometer__channel{T},
    ap2RZP::LinearMap{RotMatrix3{Float64}},
    det2RZP::LinearMap{RotMatrix3{Float64}},
)::FoV where {T <: Real}
    last_ap = ch.aperture[end]
    ap_nop = length(last_ap.outline.x1)
    ap2det = inv(det2RZP) ∘ ap2RZP
    ap = [
        ap2det(SVector(last_ap.outline.x1[ii], last_ap.outline.x2[ii], 0.0))
        for
        ii ∈ 1:ap_nop
    ]
    det_nop = length(ch.detector.outline.x1)
    det = [
        SVector(ch.detector.outline.x1[ii], ch.detector.outline.x2[ii], 0.0) for
        ii ∈ 1:det_nop
    ]

    ap_c = mean(ap)
    det_c = mean(det)
    ap_n = ap2det(SVector(0.0, 0.0, 1.0)) # get_plane_normal(ap[1], ap[2], ap[3])
    det_n = SVector(0.0, 0.0, 1.0) # get_plane_normal(det[1], det[2], det[3])
    ap_n = sign(dot(ap_n, det_n)) * ap_n
    ap_r = maximum([norm(ap_p - ap_c) for ap_p ∈ ap])
    det_r = maximum([norm(det_p - det_c) for det_p ∈ det])

    # Compute coordinate transform
    # First translate to form detector centroid as origin
    T1 = Translation(-det_c)

    # Then rotate so that aperture normal is in x-z plane
    if ap_n == SVector(0.0, 0.0, 1.0)
        R1 = LinearMap(one(RotMatrix{3, Float64}))
    else
        R1 = LinearMap(RotZ(-atan(ap_n[2] / ap_n[1])))
    end
    R1_T1 = R1 ∘ T1

    # Transform aperture centroid
    ap_tr_c = R1_T1(ap_c)
    # Rotate aperture normals
    ap_tr_n = R1(ap_n)

    # Now working on the projection to x-z plane by rotating aperture center to
    # the x-z plane but keeping the normal the same
    # This is an approximation that slighlty increases the field of view
    R_ap_c_to_xz_axis = normalize(cross(ap_tr_c, SVector(0.0, 1.0, 0.0)))
    R_ap_c_to_xz_angle = -asin(dot(ap_tr_c, SVector(0.0, 1.0, 0.0)) / norm(ap_tr_c))

    R_ap_c_to_xz = LinearMap(
        AngleAxis(
            R_ap_c_to_xz_angle,
            R_ap_c_to_xz_axis[1],
            R_ap_c_to_xz_axis[2],
            R_ap_c_to_xz_axis[3],
        ),
    )
    ap_xz_pll = SVector(ap_tr_n[3], -ap_tr_n[1])
    ap_tr_c_rotated_to_xz = R_ap_c_to_xz(ap_tr_c)
    ap_xz_c = SVector(ap_tr_c_rotated_to_xz[1], ap_tr_c_rotated_to_xz[3])
    ap_xz_ep1 = ap_xz_c + ap_r * ap_xz_pll
    ap_xz_ep2 = ap_xz_c - ap_r * ap_xz_pll
    det_xz_ep1 = SVector(-det_r, 0.0)
    det_xz_ep2 = SVector(det_r, 0.0)

    extreme_ray_11 = get_line(det_xz_ep1, ap_xz_ep1)
    extreme_ray_12 = get_line(det_xz_ep1, ap_xz_ep2)
    extreme_ray_21 = get_line(det_xz_ep2, ap_xz_ep1)
    extreme_ray_22 = get_line(det_xz_ep2, ap_xz_ep2)

    # Find the two lines that form widest field of view
    angle_11_22 = get_angle(extreme_ray_11, extreme_ray_22)
    angle_12_21 = get_angle(extreme_ray_12, extreme_ray_21)
    if angle_11_22 > angle_12_21
        fov_ray_1 = extreme_ray_11
        fov_ray_2 = extreme_ray_22
        fov_ha = angle_11_22 / 2
    else
        fov_ray_1 = extreme_ray_12
        fov_ray_2 = extreme_ray_21
        fov_ha = angle_12_21 / 2
    end

    bis_1, bis_2 = get_angle_bisector(fov_ray_1, fov_ray_2)
    fov_vertex = compute_intersection(bis_1, bis_2)
    if !isnothing(compute_intersection(ap_xz_ep1, ap_xz_ep2, bis_1))
        fov_axis_eq = bis_1
    else
        fov_axis_eq = bis_2
    end

    fov_axis = normalize(SVector(-fov_axis_eq[2], 0.0, fov_axis_eq[1]))
    fov_vertex = SVector(fov_vertex[1], 0.0, fov_vertex[2])

    # Rotate the axis and vertex back the amount aperture center was rotated
    inv_R_ap_c_to_xz = inv(R_ap_c_to_xz)
    fov_axis = inv_R_ap_c_to_xz(fov_axis) * sign(fov_axis[3])
    fov_vertex = inv_R_ap_c_to_xz(fov_vertex)

    # Add transform to rotate coordinates such that fov_axis is aligned with z-axis
    # and fov_vertex is at origin
    T2 = Translation(-fov_vertex)
    R3_axis = cross(fov_axis, SVector(0.0, 0.0, 1.0))
    R3_angle = acos(dot(fov_axis, SVector(0.0, 0.0, 1.0)))
    if R3_angle == 0 || R3_axis == SVector(0.0, 0.0, 0.0)
        R3 = LinearMap(one(RotMatrix{3, Float64}))
    else
        R3 = LinearMap(AngleAxis(R3_angle, R3_axis[1], R3_axis[2], R3_axis[3]))
    end

    # Transformation to go from detector frame to field of view frame
    det2fov = R3 ∘ T2 ∘ R1_T1

    # Transformation to go from field of view frame to R, Z, Phi frame
    # composition of fov2det (which is inv(det2fov)) and then det2RZP
    fov2RZP = det2RZP ∘ inv(det2fov)

    return FoV(fov_ha, fov2RZP)
end

"""
    get_line(
        p1::SVector{2, Float64},
        p2::SVector{2, Float64},
    )::Union{SVector{3, Float64}, Nothing}

Compute the line equation from two points in 2D space. The line is represented by the
equation ``a x + b y + c = 0`` and stored as `SVector{3, Float64}(a, b, c)`. Returns
nothing if the two points are the same.
"""
function get_line(
    p1::SVector{2, Float64},
    p2::SVector{2, Float64},
)::Union{SVector{3, Float64}, Nothing}
    dif = p2 - p1

    if dif == SVector(0.0, 0.0)
        return nothing
    end
    n = SVector(-dif[2], dif[1])
    c = -dot(p1, n)
    return SVector(n[1], n[2], c)
end

"""
    get_angle(l1::SVector{3, Float64}, l2::SVector{3, Float64})::Float64

Compute the angle in radians between two lines in 2D space. The lines are represented
by the equation ``a x + b y + c = 0`` and stored as `SVector{3, Float64}(a, b, c)`.
"""
function get_angle(l1::SVector{3, Float64}, l2::SVector{3, Float64})::Float64
    return acos(dot(l1[1:2], l2[1:2]) / (norm(l1[1:2]) * norm(l2[1:2])))
end

"""
    get_angle_bisector(
        l1::SVector{3, Float64},
        l2::SVector{3, Float64},
    )::Tuple{SVector{3, Float64}, SVector{3, Float64}}

Compute the angle bisector of two lines in 2D space. The lines are represented by the
equation ``a x + b y + c = 0`` and stored as `SVector{3, Float64}(a, b, c)`. Returns
Tuple of the anglie bisector equations in same format.
"""
function get_angle_bisector(
    l1::SVector{3, Float64},
    l2::SVector{3, Float64},
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    n1 = norm(l1[1:2])
    n2 = norm(l2[1:2])
    a1, b1, c1 = l1
    a2, b2, c2 = l2
    ab1 = SVector(a1 * n2 - a2 * n1, b1 * n2 - b2 * n1, c1 * n2 - c2 * n1)
    ab2 = SVector(a1 * n2 + a2 * n1, b1 * n2 + b2 * n1, c1 * n2 + c2 * n1)
    return ab1, ab2
end

"""
    compute_intersection(
        line1::SVector{3, Float64},
        line2::SVector{3, Float64},
    )::Union{SVector{2, Float64}, Nothing}

Compute intersection point of two infinite lines. Return nothing if they are parallel.
The lines are represented by the equation ``a x + b y + c = 0`` and stored as
`SVector{3, Float64}(a, b, c)`.
"""
function compute_intersection(
    line1::SVector{3, Float64},
    line2::SVector{3, Float64},
)::Union{SVector{2, Float64}, Nothing}
    a1, b1, c1 = line1
    a2, b2, c2 = line2
    d = a1 * b2 - a2 * b1
    if d == 0
        return nothing
    end
    x = (b1 * c2 - b2 * c1) / d
    y = (a2 * c1 - a1 * c2) / d
    return SVector(x, y)
end

"""
    compute_intersection(
        p1::SVector{2, Float64},
        p2::SVector{2, Float64},
        line::SVector{3, Float64},
    )::Union{SVector{2, Float64}, Nothing}

Compute intersection point of a finite line segment and an infinite line. Return nothing
if they are parallel or the intersection point is outside the line segment. The line
segment is represented by two points and the infinite line is represented by the
equation ``a x + b y + c = 0`` and stored as `SVector{3, Float64}(a, b, c)`.
"""
function compute_intersection(
    p1::SVector{2, Float64},
    p2::SVector{2, Float64},
    line::SVector{3, Float64},
)::Union{SVector{2, Float64}, Nothing}
    x1, y1 = p1
    x2, y2 = p2
    a, b, c = line
    d = a * (x2 - x1) + b * (y2 - y1)
    if d == 0
        return nothing
    end
    t = -(a * x1 + b * y1 + c) / d
    if t < 0 || t > 1
        return nothing
    end
    return SVector(x1 + t * (x2 - x1), y1 + t * (y2 - y1))
end