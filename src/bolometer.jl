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
frame to the (X, Y, Z) frame. In the field of view frame, the z-axis is along the
conical axis and the origin is the apex of the cone. The positive z direction is away
from the detector and towards the plasma.

Constructors:

    FoV(ha::Float64, fov2XYZ::AffineMap{T, U})
     where {T <: Rotation{3, Float64}, U <: SVector{3, Float64}}
"""
mutable struct FoV
    var"ha"::Float64 # Half angle of conical field of view in radians
    var"fov2XYZ"::AffineMap{
        T,
        U,
    } where {T <: Rotation{3, Float64}, U <: SVector{3, Float64}}
end

"""
Defining coordinate systems

(X, Y, Z) coordinate system (XYZ) is Cartesian:
X is the major radius axis for phi = 0
Y is the major radius axis for phi = pi/2
Z is the height axis.

(R, Phi, Z) coordinate system (RPZ) is Cylindrical:
R is the major radius
Phi is the toroidal angle in radians
Z is the height axis
"""
XYZ2RPZ = CylindricalFromCartesian()
RPZ2XYZ = CartesianFromCylindrical()

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
        ap2XYZ = Array{AffineMap{RotMatrix3{Float64}, SVector{3, Float64}}}(
            undef,
            length(ch.aperture),
        )

        # Add outline to apertures and detector if not present
        for (ii, ap) ∈ enumerate(ch.aperture)
            create_outline!(ap; c2o_nop=c2o_nop)
            ap2XYZ[ii] = get_transform_to_XYZ(ap)
        end
        create_outline!(ch.detector; c2o_nop=c2o_nop)
        det2XYZ = get_transform_to_XYZ(ch.detector)

        # Compute field of view for each channel
        FoVs[ch.identifier] = get_FoV(ch, ap2XYZ[end], det2XYZ)
    end
end

"""
    get_FoV(
        ch::IMASDD.bolometer_channel,
        ap2XYZ::LinearMap{RotMatrix3{Float64}},
        det2XYZ::LinearMap{RotMatrix3{Float64}},
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
    ap2XYZ::AffineMap{RotMatrix3{Float64}, SVector{3, Float64}},
    det2XYZ::AffineMap{RotMatrix3{Float64}, SVector{3, Float64}},
)::FoV where {T <: Real}
    last_ap = ch.aperture[end]
    ap_nop = length(last_ap.outline.x1)
    ap2det = inv(det2XYZ) ∘ ap2XYZ
    # Get the aperture vertices in detector frame
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

    # Compute aperture and detector centroids, normals, and radii (enclosing circle)
    ap_c = mean(ap)
    det_c = mean(det)
    ap_n = ap2det(SVector(0.0, 0.0, 1.0)) # Aperture normal is z-axis in aperture frame
    det_n = SVector(0.0, 0.0, 1.0) # Detector normal is z-axis in detector frame
    ap_r = maximum([norm(ap_p - ap_c) for ap_p ∈ ap])
    det_r = maximum([norm(det_p - det_c) for det_p ∈ det])
    ap_n = sign(dot(ap_n, det_n)) * ap_n # Ensure ap. normal points along det. normal

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

    # Now working on the projection to x-z plane by rotating aperture centre to
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

    # Rotate the axis and vertex back the amount aperture centre was rotated
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

    # Transformation to go from field of view frame to (X, Y, Z) frame
    # composition of fov2det (which is inv(det2fov)) and then det2XYZ
    fov2XYZ = det2XYZ ∘ inv(det2fov)

    return FoV(fov_ha, fov2XYZ)
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

function area_of_polygon(vertices::Vector{SVector{2, Float64}})::Float64
    n = length(vertices)
    # Shoelace formula
    return 0.5 * sum(
        [
        vertices[ii][1] * vertices[mod1(ii + 1, n)][2] -
        vertices[mod1(ii + 1, n)][1] * vertices[ii][2]
        for
        ii ∈ 1:n
    ],
    )
end

function area(
    outline::Union{
        IMASDD.bolometer__channel___detector__outline,
        IMASDD.bolometer__channel___aperture___outline,
    },
)::Float64
    return area_of_polygon(
        [SVector(outline.x1[ii], outline.x2[ii]) for ii ∈ 1:length(outline.x1)],
    )
end

function create_outline!(
    det_or_ap::Union{
        IMASDD.bolometer__channel___aperture,
        IMASDD.bolometer__channel___detector,
    };
    c2o_nop::Int64=12,
)
    if det_or_ap.geometry_type == 1
        if IMASDD.ismissing(ap, :surface)
            det_or_ap.surface = area(det_or_ap.outline)
        end
    elseif det_or_ap.geometry_type == 2
        det_or_ap.outline.x1 = det_or_ap.radius .* cos.(range(0, 2π, c2o_nop))
        det_or_ap.outline.x2 = det_or_ap.radius .* sin.(range(0, 2π, c2o_nop))
        if IMASDD.ismissing(det_or_ap, :surface)
            det_or_ap.surface = π * det_or_ap.radius^2
        end
    elseif det_or_ap.geometry_type == 3
        det_or_ap.outline.x1 = [
            -det_or_ap.x1_width / 2,
            -det_or_ap.x1_width / 2,
            det_or_ap.x1_width / 2,
            det_or_ap.x1_width / 2,
        ]
        det_or_ap.outline.x2 = [
            -det_or_ap.x2_width / 2,
            det_or_ap.x2_width / 2,
            det_or_ap.x2_width / 2,
            -det_or_ap.x2_width / 2,
        ]
        if IMASDD.ismissing(det_or_ap, :surface)
            det_or_ap.surface = det_or_ap.x1_width * det_or_ap.x2_width
        end
    end
end

function get_transform_to_XYZ(
    det_or_ap::Union{
        IMASDD.bolometer__channel___aperture,
        IMASDD.bolometer__channel___detector,
    },
)::AffineMap{RotMatrix3{Float64}, SVector{3, Float64}}
    x1 = normalize(
        SVector{3, Float64}(
            det_or_ap.x1_unit_vector.x,
            det_or_ap.x1_unit_vector.y,
            det_or_ap.x1_unit_vector.z,
        ),
    )
    x2 = normalize(
        SVector{3, Float64}(
            det_or_ap.x2_unit_vector.x,
            det_or_ap.x2_unit_vector.y,
            det_or_ap.x2_unit_vector.z,
        ),
    )
    x3 = normalize(
        SVector{3, Float64}(
            det_or_ap.x3_unit_vector.x,
            det_or_ap.x3_unit_vector.y,
            det_or_ap.x3_unit_vector.z,
        ),
    )
    origin = RPZ2XYZ(
        Cylindrical(det_or_ap.centre.r, det_or_ap.centre.phi, det_or_ap.centre.z),
    )
    return Translation(origin) ∘ LinearMap(RotMatrix3([x1; x2; x3]))
end

function add_bolometer_detector!(
    ch::IMASDD.bolometer__channel,
    centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
    x1_uv::SVector{3, Float64},
    x2_uv::SVector{3, Float64},
    x3_uv::SVector{3, Float64},
    outline::Vector{SVector{2, Float64}},
)
    add_bolometer_det_or_ap_centre!(ch.detector, centre)
    add_bolometer_det_or_ap_uv!(ch.detector, x1_uv, x2_uv, x3_uv)
    ch.detector.geometry_type = 1
    ch.detector.outline.x1 = [p[1] for p ∈ outline]
    ch.detector.outline.x2 = [p[2] for p ∈ outline]
    return ch.detector.surface = area(ch.detector.outline)
end

function add_bolometer_detector!(
    ch::IMASDD.bolometer__channel,
    centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
    x1_uv::SVector{3, Float64},
    x2_uv::SVector{3, Float64},
    x3_uv::SVector{3, Float64},
    radius::Float64,
)
    add_bolometer_det_or_ap_centre!(ch.detector, centre)
    add_bolometer_det_or_ap_uv!(ch.detector, x1_uv, x2_uv, x3_uv)
    ch.detector.geometry_type = 2
    ch.detector.radius = radius
    return ch.detector.surface = π * radius^2
end

function add_bolometer_detector!(
    ch::IMASDD.bolometer__channel,
    centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
    x1_uv::SVector{3, Float64},
    x2_uv::SVector{3, Float64},
    x3_uv::SVector{3, Float64},
    x1_width::Float64,
    x2_width::Float64,
)
    add_bolometer_det_or_ap_centre!(ch.detector, centre)
    add_bolometer_det_or_ap_uv!(ch.detector, x1_uv, x2_uv, x3_uv)
    ch.detector.geometry_type = 3
    ch.detector.x1_width = x1_width
    ch.detector.x2_width = x2_width
    return ch.detector.surface = x1_width * x2_width
end

function add_bolometer_aperture!(
    ch::IMASDD.bolometer__channel,
    centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
    x1_uv::SVector{3, Float64},
    x2_uv::SVector{3, Float64},
    x3_uv::SVector{3, Float64},
    outline::Vector{SVector{2, Float64}},
)
    resize!(ch.aperture, length(ch.aperture) + 1)
    ap = ch.aperture[end]
    add_bolometer_det_or_ap_centre!(ap, centre)
    add_bolometer_det_or_ap_uv!(ap, x1_uv, x2_uv, x3_uv)
    ap.geometry_type = 1
    ap.outline.x1 = [p[1] for p ∈ outline]
    ap.outline.x2 = [p[2] for p ∈ outline]
    return ap.surface = area(ap.outline)
end

function add_bolometer_aperture!(
    ch::IMASDD.bolometer__channel,
    centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
    x1_uv::SVector{3, Float64},
    x2_uv::SVector{3, Float64},
    x3_uv::SVector{3, Float64},
    radius::Float64,
)
    resize!(ch.aperture, length(ch.aperture) + 1)
    ap = ch.aperture[end]
    add_bolometer_det_or_ap_centre!(ap, centre)
    add_bolometer_det_or_ap_uv!(ap, x1_uv, x2_uv, x3_uv)
    ap.geometry_type = 2
    ap.radius = radius
    return ap.surface = π * radius^2
end

function add_bolometer_aperture!(
    ch::IMASDD.bolometer__channel,
    centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
    x1_uv::SVector{3, Float64},
    x2_uv::SVector{3, Float64},
    x3_uv::SVector{3, Float64},
    x1_width::Float64,
    x2_width::Float64,
)
    resize!(ch.aperture, length(ch.aperture) + 1)
    ap = ch.aperture[end]
    add_bolometer_det_or_ap_centre!(ap, centre)
    add_bolometer_det_or_ap_uv!(ap, x1_uv, x2_uv, x3_uv)
    ap.geometry_type = 3
    ap.x1_width = x1_width
    ap.x2_width = x2_width
    return ap.surface = x1_width * x2_width
end

function add_bolometer_det_or_ap_centre!(
    det_or_ap::Union{
        IMASDD.bolometer__channel___detector,
        IMASDD.bolometer__channel___aperture,
    },
    centre::Cylindrical{Float64, Float64},
)
    det_or_ap.centre.r = centre.r
    det_or_ap.centre.phi = centre.θ
    return det_or_ap.centre.z = centre.z
end

function add_bolometer_det_or_ap_centre!(
    det_or_ap::Union{
        IMASDD.bolometer__channel___detector,
        IMASDD.bolometer__channel___aperture,
    },
    centre::SVector{3, Float64},
)
    return add_bolometer_det_or_ap_centre!(det_or_ap, XYZ2RPZ(centre))
end

function add_bolometer_det_or_ap_uv!(
    det_or_ap::Union{
        IMASDD.bolometer__channel___detector,
        IMASDD.bolometer__channel___aperture,
    },
    x1_uv::SVector{3, Float64},
    x2_uv::SVector{3, Float64},
    x3_uv::SVector{3, Float64},
)
    det_or_ap.x1_unit_vector.x = x1_uv[1]
    det_or_ap.x1_unit_vector.y = x1_uv[2]
    det_or_ap.x1_unit_vector.z = x1_uv[3]
    det_or_ap.x2_unit_vector.x = x2_uv[1]
    det_or_ap.x2_unit_vector.y = x2_uv[2]
    det_or_ap.x2_unit_vector.z = x2_uv[3]
    det_or_ap.x3_unit_vector.x = x3_uv[1]
    det_or_ap.x3_unit_vector.y = x3_uv[2]
    return det_or_ap.x3_unit_vector.z = x3_uv[3]
end
