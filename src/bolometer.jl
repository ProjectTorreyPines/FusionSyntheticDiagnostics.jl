import PhysicalConstants.CODATA2018: c_0, ε_0, m_e, m_u, e
import QuadGK: quadgk, BatchIntegrand
import IMASggd: interp, get_grid_subset, get_subset_boundary, subset_do, get_TPS_mats
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

        @nospecialize(ids::IMAS.dd)=IMAS.dd();

        overwrite::Bool=false, verbose::Bool=false, kwargs...,

    )::IMAS.dd

Add bolometer to IMAS structure using a `JSON` file or Julia `Dict` and compute the
bolometer outputs. `kwargs` are passed to [`compute_bolometer!`](@ref).
"""
function add_bolometer!(
    config::Union{String, Dict{Symbol, Any}}=default_bolometer,
    @nospecialize(ids::IMAS.dd)=IMAS.dd();
    overwrite::Bool=false, verbose::Bool=false, kwargs...,
)::IMAS.dd
    add_diagnostic!(config, :bolometer, ids; overwrite=overwrite, verbose=verbose)
    compute_bolometer!(ids; kwargs...)
    return ids
end

"""
    compute_bolometer!(

        @nospecialize(ids::IMAS.dd);

        rtol::Float64=1e-3,
        n_e_gsi::Int=5,
    )

Computed the line integrated electron density from the bolometer data present in
IDS structure for all the chords. The computation is based on the edge profile data
and core profile data present in the IDS structure.
"""
function compute_bolometer!(
    @nospecialize(ids::IMAS.dd);
    c2o_nop::Int64=12,
)
    FoVs = Dict{String, FoV}()
    det2XYZs = Dict{String, AffineMap{RotMatrix3{Float64}, SVector{3, Float64}}}()
    ap2XYZs = Dict{String, Array{AffineMap{RotMatrix3{Float64}, SVector{3, Float64}}}}()
    for ch ∈ ids.bolometer.channel
        # List of transformation to go from aperture frame to X, Y, Z frame
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
        det2XYZs[ch.identifier] = det2XYZ
        ap2XYZs[ch.identifier] = ap2XYZ
    end
    return FoVs, det2XYZs, ap2XYZs
end

"""
    get_FoV(
        ch::IMAS.bolometer__channel{T},
        ap2XYZ::AffineMap{RotMatrix3{Float64}, SVector{3, Float64}},
        det2XYZ::AffineMap{RotMatrix3{Float64}, SVector{3, Float64}},
    )::FoV where {T <: Real}

Assuming aperture and detector are arrays of points described in 3D space on some

coordinate axes and that the array of points for aperture lie in a plane and the array
of points for the detector lie in a plane. This function computes the field of view
of the detector from the aperture and returns a tuple of the half angle of the field of
view and the transformation that maps from field of view coordinates to the input
coordinate frame. In the field of view coordinate frame, the field of view is a cone
with the vertex at the origin and the axis along the z-axis, and half angle that is
returned.

First, two circles are created, centered at aperture and detector and covering entire
area of the aperture and detector in detector frame. The frame is rotated about the
z-axis to get the aperture normal to lie in x-z plane. In this frame, at y=0 plane,
the aperture and detector become two line segments. The extremeties of the two line
segments are joined to create 4 extreme rays. The pair of extreme rays creating the
the largest angle between them is then chosen as the extreme rays that would form the
cone of field of view (FoV). The angle bisector is calculated to get the axis of FoV
and half-angle of FoV is calculated. These two are then returned as FoV object which
comprises of half-angle and the transformation required to go from FoV frame to XYZ
frame. In FoV frame, the conical axis is along the z axis.
"""
function get_FoV(
    ch::IMAS.bolometer__channel{T},
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

"""
    compute_intersection(
        p1::SVector{2, Float64},
        p2::SVector{2, Float64},
        line_0::SVector{3, Float64},
        line_dir::SVector{3, Float64},
    )::Union{SVector{3, Float64}, Nothing}

Compute intersection point of a finite surface and an infinite line. Return nothing
if they are parallel or the intersection point is outside the line segment.
The finite surface is represented by two points `p1` and `p2` on x-z plane (R-Z) and
created by rotating the line segment formed between `p1` and `p2` around the z-axis.
Thus, finite 2D surface in 3 dimensions is given by points:

```
x = (p1[1] * t + p2[1] * (1 - t)) * cos(phi)
y = (p1[1] * t + p2[1] * (1 - t)) * sin(phi)
z = p1[2] * t + p2[2] * (1 - t)
```

where t is in [0, 1] and phi is in [0, 2 * pi).

The infinite line is represented by an initial point `line_0` and a
direction `line_dir`. The line is given by the equation:

```
x = line_0[1] + line_dir[1] * s
y = line_0[2] + line_dir[2] * s
z = line_0[3] + line_dir[3] * s
```

where s is a real number.

Thus the intersection point is computed by solving the following equation:

```
(p1[1] * t + p2[1] * (1 - t)) * cos(phi) = line_0[1] + line_dir[1] * s
(p1[1] * t + p2[1] * (1 - t)) * sin(phi) = line_0[2] + line_dir[2] * s
p1[2] * t + p2[2] * (1 - t) = line_0[3] + line_dir[3] * s
```

for t, phi, and s. If the solution is within the range [0, 1] for t and [0, 2 * pi] for
phi, the intersection point is computed and returned. Otherwise, nothing is returned.
"""
function compute_intersection(
    p1::SVector{2, Float64},
    p2::SVector{2, Float64},
    line_0::SVector{3, Float64},
    line_dir::SVector{3, Float64},
)::Union{Array{SVector{3, Float64}}, Nothing}
    a, b = p1
    d, f = p2
    g, k, m = line_0
    h, l, n = line_dir
    # Surface:
    # x = (d + (a-d) t) cos(phi)
    # y = (d + (a-d) t) sin(phi)
    # z = f + (b-f) t
    # Line:
    # x = g + hs
    # y = k + ls
    # z = m + ns
    # Set u = d + (a-d) t
    #     v = f + (b-f) t
    # Thus for surface:
    # x = u cos(phi)
    # y = u sin(phi)
    # z = v
    # Solving for x^2 + y^2 = u^2 with line
    # u^2 = (g + hs)^2 + (k + ls)^2
    p = g^2 + k^2
    q = 2 * g * h + 2 * k * l
    r = h^2 + l^2
    # This makes the equation u^2 = p + q s + r s^2
    # Then solving for z = v = m + ns
    if n == 0                # z = m
        if b - f == 0        # Surface is parallel to x-y plane, only grazing possible
            return nothing
        else                            # Surface can intersect with the line               
            t = (m - f) / (b - f)
        end
        z = m
        u = d + (a - d) * t
        if u == 0
            return nothing
        end
        # This makes the equation r s^2 + qs + p - u^2 = 0
        ss = quadratic_roots(r, q, p - u^2)
        if isnothing(ss)
            return nothing
        end
    else
        # s = (v - m) / n
        # Thus u^2 = p + q (v - m) / n + r (v - m)^2 / n^2
        #      u^2 = r / n^2 v^2 + (q / n - 2 m r / n^2) v + p - q m / n + m^2 r / n^2
        α = r / (n^2)
        β = (q / n) - (2 * m * r / (n^2))
        γ = p - (q * m / n) + (m^2 * r / (n^2))
        # This makes the equation u^2 = α v^2 + β v + γ
        # Now we'll expand u and v in terms of t
        δ = a - d
        σ = b - f
        # This makes u = d + δ t and v = f + σ t
        # Thus u^2 = d^2 + 2 d δ t + δ^2 t^2 = α (f + σ t)^2 + β (f + σ t) + γ
        # (δ^2 - α σ^2) t^2 + (2 d δ - 2 f α σ - β σ) t + (d^2 - α f^2 - β f - γ) = 0
        A = δ^2 - α * σ^2
        B = 2 * d * δ - 2 * f * α * σ - β * σ
        C = d^2 - α * f^2 - β * f - γ
        # This makes the equation A t^2 + B t + C = 0
        # Now just solve for t if see if a solution exists
        tt = quadratic_roots(A, B, C)
        if isnothing(tt)
            return nothing
        end
        ss = Array{Float64}(undef, 0)
        for t ∈ tt
            if t < 0 || t > 1
                continue
            end
            # Now we have a valid t, so we can compute the intersection point(s)
            v = f + (b - f) * t
            s = (v - m) / n
            push!(ss, s)
        end
    end
    intersections = Array{SVector{3, Float64}}(undef, 0)
    for s ∈ ss
        x = g + h * s
        y = k + l * s
        z = m + n * s
        push!(intersections, SVector{3, Float64}(x, y, z))
    end
    return intersections
end

"""
    quadratic_roots(a::Float64, b::Float64, c::Float64)::Union{Array{Float64}, Nothing}

Compute the roots of a quadratic equation `a x^2 + b x + c = 0`. Return nothing if the
discriminant is negative, return an array with one element if the discriminant is zero,
and return an array with two elements if the discriminant is positive.
"""
function quadratic_roots(
    a::Float64,
    b::Float64,
    c::Float64,
)::Union{Array{Float64}, Nothing}
    if a == 0
        if b == 0
            return nothing
        else
            x = -c / b
            return x
        end
    end
    delta = b^2 - 4 * a * c
    xx = Array{Float64}(undef, 0)
    if delta < 0
        return nothing
    elseif delta == 0
        x = -b / (2 * a)
        return Array{Float64}([x])
    end
    x1 = (-b + sqrt(delta)) / (2 * a)
    x2 = (-b - sqrt(delta)) / (2 * a)
    return Array{Float64}([x1, x2])
end

"""
    area_of_polygon(vertices::Vector{SVector{2, Float64}})::Float64

Function to calculate area of arbitrary non-self intersecting polygon.
"""
function area_of_polygon(vertices::Vector{SVector{2, Float64}})::Float64
    n = length(vertices)
    # Shoelace formula
    return 0.5 * abs(
        sum(
            [
            vertices[ii][1] * vertices[mod1(ii + 1, n)][2] -
            vertices[mod1(ii + 1, n)][1] * vertices[ii][2]
            for
            ii ∈ 1:n
        ],
        ),
    )
end

"""
    area(
    outline::Union{
        IMAS.bolometer__channel___detector__outline,
        IMAS.bolometer__channel___aperture___outline,
    },

)::Float64

Function to return area of an aperture or detector surface.
"""
function area(
    outline::Union{
        IMAS.bolometer__channel___detector__outline,
        IMAS.bolometer__channel___aperture___outline,
    },
)::Float64
    return area_of_polygon(
        [SVector(outline.x1[ii], outline.x2[ii]) for ii ∈ eachindex(outline.x1)],
    )
end

"""
    project_3D_line_to_z0(
        p1::SVector{3, Float64},
        p2::SVector{3, Float64},
    )::SVector{2, Float64}

Projects a 3D line between points `p1` and `p2` onto z=0 plane.
"""
function project_3D_line_to_z0(
    p1::SVector{3, Float64},
    p2::SVector{3, Float64},
)::SVector{2, Float64}
    t = -p1[3] / (p2[3] - p1[3])
    proj = p1 + t .* (p2 - p1)
    return SVector(proj[1], proj[2])
end

"""
    get_lit_sa(
        source::SVector{3, Float64},
        ch::IMAS.bolometer__channel,
        XYZ2det::AffineMap{RotMatrix3{Float64}, SVector{3, Float64}},
        ap2XYZs::Array{AffineMap{RotMatrix3{Float64}, SVector{3, Float64}}},
    )::Float64

Get illuminated solid angle from source point (in XYZ frame) on the detector through
the array of apertures of a channel. The solid angle is calculated as ratio of
illuminated area on the detector to square of distance from source to detector. The
illuminated area is calculated by clipping the detector outline with the outlines that
apertures make on the detector plane when projected from the source point.
"""
function get_lit_sa(
    source::SVector{3, Float64},
    ch::IMAS.bolometer__channel,
    XYZ2det::AffineMap{RotMatrix3{Float64}, SVector{3, Float64}},
    ap2XYZs::Array{AffineMap{RotMatrix3{Float64}, SVector{3, Float64}}},
)::Float64
    # src is in detector frame of reference
    src = XYZ2det(source)
    det_list = [
        SVector(ch.detector.outline.x1[ii], ch.detector.outline.x2[ii]) for
        ii ∈ eachindex(ch.detector.outline.x1)
    ]
    lit_region = deepcopy(det_list)
    for (ap_ind, ap) ∈ enumerate(ch.aperture)
        # ap2det is transformation to convert from aperture frame of detector frame
        ap2det = XYZ2det ∘ ap2XYZs[ap_ind]
        ap_list = [
            ap2det(SVector(ap.outline.x1[ii], ap.outline.x2[ii], 0.0)) for
            ii ∈ eachindex(ap.outline.x1)
        ]
        # clip_list stores projections of src to aperture vertex onto detector plane
        clip_list = [project_3D_line_to_z0(src, app) for app ∈ ap_list]
        # lit_region is iteratively clipped for each aperture
        lit_region = clip(lit_region, clip_list)
    end
    # In detector frame, detector center is at origin
    return area_of_polygon(lit_region) / norm(src)^2
end

"""
    clip(
        subject_list::Vector{SVector{2, Float64}},
        clip_list::Vector{SVector{2, Float64}},
    )::Vector{SVector{2, Float64}}

Clip a list of vertices that form a polygon by a list of vertices that form a convex
clipping ploygon. It is assumed the lists are ordered.
"""
function clip(
    subject_list::Vector{SVector{2, Float64}},
    clip_list::Vector{SVector{2, Float64}},
)::Vector{SVector{2, Float64}}
    N = length(clip_list)
    clip_edges = [
        get_line(clip_list[ii], clip_list[mod1(ii + 1, N)]) for
        ii ∈ eachindex(clip_list)
    ]
    for ii ∈ eachindex(clip_edges)
        if !right_of_edge(clip_list[mod1(ii + 2, N)], clip_edges[ii])
            clip_edges[ii] = -1.0 .* clip_edges[ii]
        end
    end
    return clip(subject_list, clip_edges)
end

"""
    clip(
        subject_list::Vector{SVector{2, Float64}},
        clip_list::Vector{SVector{3, Float64}},
    )::Vector{SVector{2, Float64}}

Clip a list of vertices that form a polygon by a list of edges that form a convex
clipping ploygon. It is assumed the lists are ordered and the clip_list edges are such
that all vertices of clipping polygon are on the right side of the edges.
"""
function clip(
    subject_list::Vector{SVector{2, Float64}},
    clip_list::Vector{SVector{3, Float64}},
)::Vector{SVector{2, Float64}}
    out_list = deepcopy(subject_list)
    for clip_edge ∈ clip_list
        inp_list = deepcopy(out_list)
        N = length(inp_list)
        out_list = SVector{2, Float64}[]
        for ii ∈ eachindex(inp_list)
            current_point = inp_list[ii]
            prev_point = inp_list[mod1(ii - 1, N)]
            # Compute intersection of finite line segment with infinite clip edge
            int_point = compute_intersection(prev_point, current_point, clip_edge)
            if !isnothing(int_point)
                push!(out_list, int_point)
            end
            if right_of_edge(current_point, clip_edge) > 0
                push!(out_list, current_point)
            end
        end
    end
    return out_list
end

"""
    right_of_edge(point::SVector{2,Float64}, edge::SVector{3,Float64})::Bool

Check if a point is on the right side of an edge.
"""
function right_of_edge(point::SVector{2, Float64}, edge::SVector{3, Float64})::Bool
    x, y = point
    a, b, c = edge
    return a * x + b * y + c > 0
end

"""
    create_outline!(
        det_or_ap::Union{
            IMAS.bolometer__channel___aperture,
            IMAS.bolometer__channel___detector,
        };
        c2o_nop::Int64=12,
    )

For a deetector or aperture, regardless of how the surface is defined (outline,
circular, or rectangular), this function converts the surface geometry type to 1
(outline) by computing representative outline points. For rectangular geometry,
the four corners are computed ot create outline. For circular geomtry, `c2o_nop` number
of points are used to approximate the circle with an outline of `c2o_nop` edges.

This function is used to keep all aperture and detector surfaces in outline format for
single format computing later.
"""
function create_outline!(
    det_or_ap::Union{
        IMAS.bolometer__channel___aperture,
        IMAS.bolometer__channel___detector,
    };
    c2o_nop::Int64=12,
)
    if det_or_ap.geometry_type == 1
        if IMAS.ismissing(ap, :surface)
            det_or_ap.surface = area(det_or_ap.outline)
        end
    elseif det_or_ap.geometry_type == 2
        angles = range(0, 2π; length=c2o_nop)[1:end-1]
        det_or_ap.outline.x1 = det_or_ap.radius .* cos.(angles)
        det_or_ap.outline.x2 = det_or_ap.radius .* sin.(angles)
        if IMAS.ismissing(det_or_ap, :surface)
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
        if IMAS.ismissing(det_or_ap, :surface)
            det_or_ap.surface = det_or_ap.x1_width * det_or_ap.x2_width
        end
    end
end

"""
    get_transform_to_XYZ(
        det_or_ap::Union{
            IMAS.bolometer__channel___aperture,
            IMAS.bolometer__channel___detector,
        },
    )::AffineMap{RotMatrix3{Float64}, SVector{3, Float64}}

For an aperture(or a detctor), use the stored unit vectors to define to transform to
go from aperture (or detector) frame to machine's XYZ frame.

First rotate to become parallel to XYZ frame and then translate by the XYZ position of
origin.
"""
function get_transform_to_XYZ(
    det_or_ap::Union{
        IMAS.bolometer__channel___aperture,
        IMAS.bolometer__channel___detector,
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

"""
    add_bolometer_detector!(
        ch::IMAS.bolometer__channel,
        centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
        x1_uv::SVector{3, Float64},
        x2_uv::SVector{3, Float64},
        x3_uv::SVector{3, Float64},
        outline::Vector{SVector{2, Float64}},
    )

Add an arbitrary detector to a bolometer channel with supplied unit vectors and outline
provided as an array of points.
"""
function add_bolometer_detector!(
    ch::IMAS.bolometer__channel,
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

"""
    add_bolometer_detector!(
        ch::IMAS.bolometer__channel,
        centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
        x1_uv::SVector{3, Float64},
        x2_uv::SVector{3, Float64},
        x3_uv::SVector{3, Float64},
        radius::Float64,
    )

Add a circular detector to a bolometer channel with supplied unit vectors.
"""
function add_bolometer_detector!(
    ch::IMAS.bolometer__channel,
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

"""
    add_bolometer_detector!(
        ch::IMAS.bolometer__channel,
        centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
        x1_uv::SVector{3, Float64},
        x2_uv::SVector{3, Float64},
        x3_uv::SVector{3, Float64},
        x1_width::Float64,
        x2_width::Float64,
    )

Add a rectangular detector to a bolometer channel with supplied unit vectors.
"""
function add_bolometer_detector!(
    ch::IMAS.bolometer__channel,
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

"""
    add_bolometer_aperture!(
        ch::IMAS.bolometer__channel,
        centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
        x1_uv::SVector{3, Float64},
        x2_uv::SVector{3, Float64},
        x3_uv::SVector{3, Float64},
        outline::Vector{SVector{2, Float64}},
    )

Add an arbitrary aperture to a bolometer channel with supplied unit vectors and outline
provided as an array of points.
"""
function add_bolometer_aperture!(
    ch::IMAS.bolometer__channel,
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

"""
    add_bolometer_aperture!(
        ch::IMAS.bolometer__channel,
        centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
        x1_uv::SVector{3, Float64},
        x2_uv::SVector{3, Float64},
        x3_uv::SVector{3, Float64},
        radius::Float64,
    )

Add a circular aperture to a bolometer channel with supplied unit vectors.
"""
function add_bolometer_aperture!(
    ch::IMAS.bolometer__channel,
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

"""
    add_bolometer_aperture!(
        ch::IMAS.bolometer__channel,
        centre::Union{Cylindrical{Float64, Float64}, SVector{3, Float64}},
        x1_uv::SVector{3, Float64},
        x2_uv::SVector{3, Float64},
        x3_uv::SVector{3, Float64},
        x1_width::Float64,
        x2_width::Float64,
    )

Add a rectangular aperture to a bolometer channel with supplied unit vectors.
"""
function add_bolometer_aperture!(
    ch::IMAS.bolometer__channel,
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

"""
    add_bolometer_det_or_ap_centre!(
        det_or_ap::Union{
            IMAS.bolometer__channel___detector,
            IMAS.bolometer__channel___aperture,
        },
        centre::Cylindrical{Float64, Float64},
    )

Add detector or aperture center to bolometer__channel___aperture or
bolometer__channel___detector if center provided in Cartesian coordinates.
"""
function add_bolometer_det_or_ap_centre!(
    det_or_ap::Union{
        IMAS.bolometer__channel___detector,
        IMAS.bolometer__channel___aperture,
    },
    centre::Cylindrical{Float64, Float64},
)
    det_or_ap.centre.r = centre.r
    det_or_ap.centre.phi = centre.θ
    return det_or_ap.centre.z = centre.z
end

"""
    add_bolometer_det_or_ap_centre!(
        det_or_ap::Union{
            IMAS.bolometer__channel___detector,
            IMAS.bolometer__channel___aperture,
        },
        centre::SVector{3, Float64},
    )

Add detector or aperture center to bolometer__channel___aperture or
bolometer__channel___detector if center provided in Cartesian coordinates.
"""
function add_bolometer_det_or_ap_centre!(
    det_or_ap::Union{
        IMAS.bolometer__channel___detector,
        IMAS.bolometer__channel___aperture,
    },
    centre::SVector{3, Float64},
)
    return add_bolometer_det_or_ap_centre!(det_or_ap, XYZ2RPZ(centre))
end

"""
    add_bolometer_det_or_ap_uv!(
        det_or_ap::Union{
            IMAS.bolometer__channel___detector,
            IMAS.bolometer__channel___aperture,
        },
        x1_uv::SVector{3, Float64},
        x2_uv::SVector{3, Float64},
        x3_uv::SVector{3, Float64},
    )

Add unit vectors to bolometer__channel___aperture or bolometer__channel___detector
"""
function add_bolometer_det_or_ap_uv!(
    det_or_ap::Union{
        IMAS.bolometer__channel___detector,
        IMAS.bolometer__channel___aperture,
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
