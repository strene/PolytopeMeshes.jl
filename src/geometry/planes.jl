# Planes, rotations and convex hulls in 3D.

"""
    normal_from_points(p)

Unit normal of the plane through the (approximately coplanar) points `p`.
"""
function normal_from_points(p::AbstractVector{Vec3})
    length(p) >= 3 || throw(ArgumentError("At least three points are needed to define a plane"))
    c = sum(p) / length(p)
    n = zero(Vec3)
    k = length(p)
    for i in 1:k
        n += cross(p[i] - c, p[mod1(i + 1, k)] - c)
    end
    if norm(n) < 1e-14
        # Points not ordered as a polygon: use the best pair of vectors
        best = 0.0
        for i in 2:k, j in (i + 1):k
            v = cross(p[i] - p[1], p[j] - p[1])
            if norm(v) > best
                best = norm(v)
                n = v
            end
        end
    end
    norm(n) > 0 || throw(ArgumentError("Points are collinear"))
    return n / norm(n)
end

"""
    plane_frame(p)

Return `(R, c)` where `c` is the mean of the planar points `p` and `R` is an
orthogonal matrix with rows `(u, v, n)` such that `R * (x - c)` maps the plane
to `z = 0`.
"""
function plane_frame(p::AbstractVector{Vec3})
    n = normal_from_points(p)
    c = sum(p) / length(p)
    # Pick the first in-plane axis along the first polygon edge if possible
    u = p[2] - p[1]
    u = u - dot(u, n) * n
    if norm(u) < 1e-12
        u = abs(n[1]) < 0.9 ? cross(n, Vec3(1, 0, 0)) : cross(n, Vec3(0, 1, 0))
    end
    u = u / norm(u)
    v = cross(n, u)
    R = SMatrix{3, 3, Float64}(u[1], v[1], n[1], u[2], v[2], n[2], u[3], v[3], n[3])
    return R, c
end

to_plane(R, c, x::Vec3) = (y = R * (x - c); Vec2(y[1], y[2]))
from_plane(R, c, x::Vec2) = R' * Vec3(x[1], x[2], 0.0) + c

"""
    rotation_matrix(theta, axis)

Rotation by angle `theta` around `axis` (Rodrigues' formula).
"""
function rotation_matrix(theta, axis)
    a = Vec3(axis) / norm(axis)
    K = SMatrix{3, 3, Float64}(0, a[3], -a[2], -a[3], 0, a[1], a[2], -a[1], 0)
    return SMatrix{3, 3, Float64}(I) + sin(theta) * K + (1 - cos(theta)) * K * K
end

"""
    convex_hull_planes(points; tol)

Facet planes of the convex hull of `points` as vectors of unit outward normals
and offsets, such that the hull is `{x : dot(n, x) <= d}`. Uses a brute force
algorithm meant for a small number of points.
"""
function convex_hull_planes(points::AbstractVector{Vec3}; tol = 1e-10)
    lo, hi = bounding_box(points)
    L = norm(hi - lo)
    np = length(points)
    normals = Vec3[]
    offsets = Float64[]
    for i in 1:np, j in (i + 1):np, k in (j + 1):np
        n = cross(points[j] - points[i], points[k] - points[i])
        nn = norm(n)
        nn < 1e-12 * L^2 && continue
        n = n / nn
        d = dot(n, points[i])
        s = [dot(n, p) - d for p in points]
        if all(s .<= tol * L)
            # n is outward
        elseif all(s .>= -tol * L)
            n = -n
            d = -d
        else
            continue
        end
        dup = false
        for (m, e) in zip(normals, offsets)
            if dot(m, n) > 1 - 1e-10 && abs(e - d) <= tol * L
                dup = true
                break
            end
        end
        if !dup
            push!(normals, n)
            push!(offsets, d)
        end
    end
    return normals, offsets
end
